#!/usr/bin/env python
"""
Run anomaly detection on protein and DNA embeddings to identify novel pathogens.
Supports multiple anomaly detection methods and dimensionality reduction techniques.
Added support for Isolation Forest, VAE, and DeepSVDD integrated models.
"""

import argparse
import numpy as np
import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # For headless environments
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap
import logging
import warnings
import json
from joblib import dump, load, Parallel, delayed
from collections import Counter
import torch
warnings.filterwarnings("ignore")

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('anomaly_detection')

def parse_args():
    parser = argparse.ArgumentParser(description='Run anomaly detection on embeddings')
    
    # Input files
    parser.add_argument('--protein_embeddings', help='Protein embeddings NPY file')
    parser.add_argument('--protein_ids', help='Protein IDs text file')
    parser.add_argument('--protein_annotations', help='Protein taxonomy annotations')
    parser.add_argument('--protein_anomalies', help='Pre-computed protein anomalies JSON file (from embedding step)')
    
    parser.add_argument('--dna_embeddings', help='DNA embeddings NPY file')
    parser.add_argument('--contig_ids', help='Contig IDs text file')
    parser.add_argument('--contig_annotations', help='Contig taxonomy annotations')
    parser.add_argument('--dna_anomalies', help='Pre-computed DNA anomalies JSON file (from embedding step)')
    
    # Output files
    parser.add_argument('--output_protein', help='Output file for protein anomaly results')
    parser.add_argument('--output_dna', help='Output file for DNA anomaly results')
    parser.add_argument('--save_models', help='Directory to save trained models')
    parser.add_argument('--visualize', help='Directory to save visualization plots')
    
    # Method selection
    parser.add_argument('--methods', default='isolation_forest,vae,dbscan',
                      help='Comma-separated list of anomaly detection methods')
    
    # Parameters
    parser.add_argument('--contamination', type=float, default=0.05,
                      help='Expected fraction of anomalies')
    parser.add_argument('--dimensions_reduction', default='umap,tsne',
                      help='Comma-separated list of dimension reduction methods')
    parser.add_argument('--perplexity', type=int, default=30,
                      help='Perplexity parameter for t-SNE')
    parser.add_argument('--random_state', type=int, default=42, 
                      help='Random state for reproducibility')
    parser.add_argument('--threads', type=int, default=8,
                      help='Number of CPU threads')
    
    # VAE/DeepSVDD specific parameters
    parser.add_argument('--latent_dim', type=int, default=64,
                      help='Latent dimension size for VAE or DeepSVDD')
    parser.add_argument('--epochs', type=int, default=100,
                      help='Number of training epochs for neural network models')
    parser.add_argument('--learning_rate', type=float, default=0.001,
                      help='Learning rate for neural network training')
    parser.add_argument('--load_models', action='store_true',
                      help='Load pre-trained models if available')
    
    # Ensemble parameters
    parser.add_argument('--ensemble', action='store_true',
                      help='Use ensemble of multiple methods')
    parser.add_argument('--ensemble_threshold', type=float, default=0.5,
                      help='Threshold for ensemble voting (fraction of methods)')
    
    return parser.parse_args()

def load_data(embedding_file, ids_file, annotation_file=None, anomalies_file=None):
    """Load embeddings, IDs, and optional annotations and pre-computed anomalies."""
    # Load embeddings
    embeddings = np.load(embedding_file)
    
    # Load IDs
    with open(ids_file, 'r') as f:
        ids = [line.strip() for line in f]
    
    # Load annotations if provided
    annotations = None
    if annotation_file:
        try:
            anno_df = pd.read_csv(annotation_file, sep='\t')
            # Create a dictionary from the annotation file
            if 'protein_id' in anno_df.columns:
                id_col = 'protein_id'
            elif 'contig_id' in anno_df.columns:
                id_col = 'contig_id'
            else:
                id_col = anno_df.columns[0]
            
            annotations = dict(zip(anno_df[id_col], anno_df['taxonomy']))
        except Exception as e:
            logger.warning(f"Error loading annotations: {e}")
            annotations = None
    
    # Load pre-computed anomalies if provided
    precomputed_anomalies = None
    if anomalies_file and os.path.exists(anomalies_file):
        try:
            with open(anomalies_file, 'r') as f:
                anomalies_data = json.load(f)
            
            # Extract anomaly labels and scores
            anomaly_labels = {}
            anomaly_scores = {}
            
            for result in anomalies_data.get('results', []):
                id = result.get('id')
                if id:
                    anomaly_labels[id] = result.get('is_anomaly', False)
                    anomaly_scores[id] = result.get('score', 0.0)
            
            precomputed_anomalies = {
                'method': anomalies_data.get('method', 'unknown'),
                'contamination': anomalies_data.get('contamination', 0.05),
                'labels': anomaly_labels,
                'scores': anomaly_scores
            }
            
            logger.info(f"Loaded pre-computed anomalies: {sum(anomaly_labels.values())} anomalies found")
        except Exception as e:
            logger.warning(f"Error loading pre-computed anomalies: {e}")
    
    return embeddings, ids, annotations, precomputed_anomalies

def reduce_dimensions(embeddings, methods, perplexity=30, n_components=2, random_state=42, n_jobs=8):
    """Apply dimensionality reduction to embeddings."""
    results = {}
    
    # Apply standardization
    logger.info("Standardizing embeddings")
    scaler = StandardScaler()
    scaled_embeddings = scaler.fit_transform(embeddings)
    
    # Apply PCA as initial dimension reduction for high-dimensional embeddings
    if embeddings.shape[1] > 100:
        logger.info("Applying PCA pre-processing to reduce dimensions to 50")
        pca_50 = PCA(n_components=50, random_state=random_state)
        embeddings_reduced = pca_50.fit_transform(scaled_embeddings)
    else:
        embeddings_reduced = scaled_embeddings
    
    # Apply requested dimension reduction methods
    for method in methods:
        logger.info(f"Applying {method} dimensionality reduction")
        
        if method.lower() == 'pca':
            reducer = PCA(n_components=n_components, random_state=random_state)
            reduced = reducer.fit_transform(scaled_embeddings)
            results['pca'] = reduced
            
        elif method.lower() == 'tsne':
            tsne = TSNE(n_components=n_components, perplexity=perplexity, 
                         n_jobs=n_jobs, random_state=random_state)
            reduced = tsne.fit_transform(embeddings_reduced)  # Use the preprocessed embeddings
            results['tsne'] = reduced
            
        elif method.lower() == 'umap':
            reducer = umap.UMAP(n_components=n_components, random_state=random_state,
                                n_neighbors=min(15, len(embeddings)-1))
            reduced = reducer.fit_transform(embeddings_reduced)  # Use the preprocessed embeddings
            results['umap'] = reduced
    
    return results, scaler

def isolation_forest_detection(embeddings, contamination=0.05, random_state=42, n_jobs=8):
    """Detect anomalies using Isolation Forest."""
    from sklearn.ensemble import IsolationForest
    
    logger.info("Running Isolation Forest anomaly detection")
    
    # Initialize and fit the model
    iso_forest = IsolationForest(
        contamination=contamination,
        random_state=random_state,
        n_jobs=n_jobs,
        verbose=0
    )
    
    # Fit and predict
    y_pred = iso_forest.fit_predict(embeddings)
    
    # Convert to binary labels and anomaly scores
    # In Isolation Forest: -1 for anomalies, 1 for normal
    # We convert to: 1 for anomalies, 0 for normal
    labels = np.where(y_pred == -1, 1, 0)
    scores = -iso_forest.decision_function(embeddings)  # Higher score means more anomalous
    
    return labels, scores, iso_forest

def dbscan_detection(embeddings, contamination=0.05):
    """Detect anomalies using DBSCAN."""
    from sklearn.cluster import DBSCAN
    from sklearn.neighbors import NearestNeighbors
    
    logger.info("Running DBSCAN anomaly detection")
    
    # Determine epsilon parameter using nearest neighbors distances
    k = max(5, int(len(embeddings) * 0.05))  # At least 5 neighbors, or 5% of data
    k = min(k, len(embeddings) - 1)  # Can't have more neighbors than points - 1
    
    nbrs = NearestNeighbors(n_neighbors=k).fit(embeddings)
    distances, _ = nbrs.kneighbors(embeddings)
    
    # Sort distances to k-th neighbor
    dist_desc = sorted(distances[:, -1], reverse=True)
    
    # Choose epsilon as the distance where the slope changes most dramatically
    # or use a percentile if the data is too small
    if len(embeddings) > 100:
        slopes = [dist_desc[i] - dist_desc[i+1] for i in range(len(dist_desc)-1)]
        eps = dist_desc[slopes.index(max(slopes))]
    else:
        eps = np.percentile(dist_desc, 75)  # 75th percentile
    
    # Run DBSCAN
    clustering = DBSCAN(eps=eps, min_samples=k).fit(embeddings)
    
    # Get labels (-1 is noise/anomaly in DBSCAN)
    labels = np.where(clustering.labels_ == -1, 1, 0)
    
    # Compute anomaly scores using distance to nearest cluster
    if np.sum(labels) == len(embeddings):  # All points are anomalies
        scores = np.ones(len(embeddings))
    else:
        # Compute distance to nearest non-anomaly point for each point
        normal_points = embeddings[labels == 0]
        nbrs = NearestNeighbors(n_neighbors=1).fit(normal_points)
        distances, _ = nbrs.kneighbors(embeddings)
        scores = distances.flatten()
    
    # Adjust anomaly count if it's too different from the expected contamination
    if abs(np.mean(labels) - contamination) > 0.1:
        # Too few anomalies - mark the most anomalous points
        if np.mean(labels) < contamination:
            threshold = np.percentile(scores, 100 * (1 - contamination))
            labels = np.where(scores >= threshold, 1, 0)
        # Too many anomalies - keep only the most anomalous ones
        else:
            anomaly_indices = np.where(labels == 1)[0]
            anomaly_scores = scores[anomaly_indices]
            
            # Sort anomalies by score and keep top contamination %
            n_to_keep = int(contamination * len(embeddings))
            threshold_idx = np.argsort(anomaly_scores)[-n_to_keep:]
            
            new_labels = np.zeros_like(labels)
            new_labels[anomaly_indices[threshold_idx]] = 1
            labels = new_labels
    
    return labels, scores, None  # No model to return for DBSCAN

def local_outlier_factor_detection(embeddings, contamination=0.05, n_jobs=8):
    """Detect anomalies using Local Outlier Factor."""
    from sklearn.neighbors import LocalOutlierFactor
    
    logger.info("Running Local Outlier Factor anomaly detection")
    
    # Initialize and fit the model
    lof = LocalOutlierFactor(
        n_neighbors=20,
        contamination=contamination,
        n_jobs=n_jobs
    )
    
    # Fit and predict
    y_pred = lof.fit_predict(embeddings)
    
    # Convert to binary labels
    # In LOF: -1 for anomalies, 1 for normal
    # We convert to: 1 for anomalies, 0 for normal
    labels = np.where(y_pred == -1, 1, 0)
    
    # Negative outlier factor (higher = more anomalous)
    scores = -lof.negative_outlier_factor_
    
    return labels, scores, None  # LOF doesn't provide a model that can be saved

def one_class_svm_detection(embeddings, contamination=0.05):
    """Detect anomalies using One-Class SVM."""
    from sklearn.svm import OneClassSVM
    
    logger.info("Running One-Class SVM anomaly detection")
    
    # Initialize and fit the model
    ocsvm = OneClassSVM(nu=contamination, kernel='rbf', gamma='scale')
    
    # Fit and predict
    y_pred = ocsvm.fit_predict(embeddings)
    
    # Convert to binary labels
    # In OneClassSVM: -1 for anomalies, 1 for normal
    # We convert to: 1 for anomalies, 0 for normal
    labels = np.where(y_pred == -1, 1, 0)
    
    # Decision function (lower = more anomalous)
    scores = -ocsvm.decision_function(embeddings)
    
    return labels, scores, ocsvm

def setup_vae_model(input_dim, latent_dim, device):
    """Set up VAE model for anomaly detection."""
    try:
        import torch.nn as nn
        import torch.nn.functional as F
        
        class VAE(nn.Module):
            def __init__(self, input_dim, hidden_dim=256, latent_dim=64):
                super(VAE, self).__init__()
                
                # Encoder
                self.fc1 = nn.Linear(input_dim, hidden_dim)
                self.bn1 = nn.BatchNorm1d(hidden_dim)
                self.fc2 = nn.Linear(hidden_dim, hidden_dim)
                self.bn2 = nn.BatchNorm1d(hidden_dim)
                self.fc31 = nn.Linear(hidden_dim, latent_dim)  # mean
                self.fc32 = nn.Linear(hidden_dim, latent_dim)  # log variance
                
                # Decoder
                self.fc4 = nn.Linear(latent_dim, hidden_dim)
                self.bn3 = nn.BatchNorm1d(hidden_dim)
                self.fc5 = nn.Linear(hidden_dim, hidden_dim)
                self.bn4 = nn.BatchNorm1d(hidden_dim)
                self.fc6 = nn.Linear(hidden_dim, input_dim)
                
                self.latent_dim = latent_dim
            
            def encode(self, x):
                h1 = F.relu(self.bn1(self.fc1(x)))
                h2 = F.relu(self.bn2(self.fc2(h1)))
                return self.fc31(h2), self.fc32(h2)
            
            def reparameterize(self, mu, logvar):
                std = torch.exp(0.5 * logvar)
                eps = torch.randn_like(std)
                return mu + eps * std
            
            def decode(self, z):
                h4 = F.relu(self.bn3(self.fc4(z)))
                h5 = F.relu(self.bn4(self.fc5(h4)))
                return self.fc6(h5)
            
            def forward(self, x):
                mu, logvar = self.encode(x)
                z = self.reparameterize(mu, logvar)
                return self.decode(z), mu, logvar
            
        logger.info(f"Creating VAE model with input dim {input_dim} and latent dim {latent_dim}")
        model = VAE(input_dim, hidden_dim=256, latent_dim=latent_dim)
        model = model.to(device)
        return model
    
    except ImportError:
        logger.error("PyTorch not properly installed. Check your installation.")
        sys.exit(1)

def setup_deep_svdd_model(input_dim, latent_dim, device):
    """Set up Deep SVDD model for anomaly detection."""
    try:
        import torch.nn as nn
        import torch.nn.functional as F
        
        class DeepSVDD(nn.Module):
            def __init__(self, input_dim, hidden_dim=256, output_dim=64):
                super(DeepSVDD, self).__init__()
                self.fc1 = nn.Linear(input_dim, hidden_dim)
                self.bn1 = nn.BatchNorm1d(hidden_dim)
                self.fc2 = nn.Linear(hidden_dim, hidden_dim)
                self.bn2 = nn.BatchNorm1d(hidden_dim)
                self.fc3 = nn.Linear(hidden_dim, output_dim)
                
                # Initialize center to None (will be set during training)
                self.center = None
                
            def forward(self, x):
                x = F.relu(self.bn1(self.fc1(x)))
                x = F.relu(self.bn2(self.fc2(x)))
                x = self.fc3(x)
                return x
            
            def init_center(self, embeddings, eps=0.1):
                """Initialize hypersphere center as the mean of all embeddings."""
                self.center = torch.mean(embeddings, dim=0)
                # Add small constant to avoid collapse to a single point
                self.center[(abs(self.center) < eps) & (self.center < 0)] = -eps
                self.center[(abs(self.center) < eps) & (self.center >= 0)] = eps
        
        logger.info(f"Creating DeepSVDD model with input dim {input_dim} and output dim {latent_dim}")
        model = DeepSVDD(input_dim, hidden_dim=256, output_dim=latent_dim)
        model = model.to(device)
        return model
    
    except ImportError:
        logger.error("PyTorch not properly installed. Check your installation.")
        sys.exit(1)

def train_vae_model(model, dataloader, epochs=100, lr=0.001, device='cpu'):
    """Train VAE model."""
    import torch
    import torch.optim as optim
    import torch.nn.functional as F
    
    optimizer = optim.Adam(model.parameters(), lr=lr)
    model.train()
    
    logger.info(f"Training VAE model for {epochs} epochs")
    for epoch in range(epochs):
        total_loss = 0
        
        for batch_x in dataloader:
            batch_x = batch_x.to(device)
            optimizer.zero_grad()
            
            # Forward pass
            recon_x, mu, logvar = model(batch_x)
            
            # Compute loss
            recon_loss = F.mse_loss(recon_x, batch_x, reduction='sum')
            kl_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
            loss = recon_loss + kl_loss
            
            # Backward pass
            loss.backward()
            optimizer.step()
            
            total_loss += loss.item()
        
        # Log progress every 10 epochs
        if (epoch + 1) % 10 == 0:
            logger.info(f"Epoch {epoch+1}/{epochs}, Loss: {total_loss / len(dataloader.dataset):.4f}")
    
    model.eval()
    return model

def train_deep_svdd_model(model, dataloader, epochs=100, lr=0.001, device='cpu'):
    """Train Deep SVDD model."""
    import torch
    import torch.optim as optim
    
    # First pass to initialize center
    logger.info("Initializing DeepSVDD center")
    model.eval()
    embeddings_list = []
    with torch.no_grad():
        for batch_x in dataloader:
            batch_x = batch_x.to(device)
            embeddings = model(batch_x)
            embeddings_list.append(embeddings)
    
    # Concatenate all embeddings and initialize center
    all_embeddings = torch.cat(embeddings_list, dim=0)
    model.init_center(all_embeddings)
    
    # Train model
    optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=1e-6)
    model.train()
    
    logger.info(f"Training DeepSVDD model for {epochs} epochs")
    for epoch in range(epochs):
        total_loss = 0
        
        for batch_x in dataloader:
            batch_x = batch_x.to(device)
            optimizer.zero_grad()
            
            # Forward pass
            embeddings = model(batch_x)
            
            # Compute distance to center (loss)
            dist = torch.sum((embeddings - model.center) ** 2, dim=1)
            loss = torch.mean(dist)
            
            # Backward pass
            loss.backward()
            optimizer.step()
            
            total_loss += loss.item() * batch_x.size(0)
        
        # Log progress every 10 epochs
        if (epoch + 1) % 10 == 0:
            logger.info(f"Epoch {epoch+1}/{epochs}, Loss: {total_loss / len(dataloader.dataset):.4f}")
    
    model.eval()
    return model

def variational_autoencoder_detection(embeddings, contamination=0.05, random_state=42, 
                                     latent_dim=64, epochs=100, lr=0.001, batch_size=32, device='cpu',
                                     load_model_path=None):
    """Detect anomalies using a Variational Autoencoder."""
    try:
        import torch
        import torch.utils.data as data_utils
        
        # Set random seeds for reproducibility
        torch.manual_seed(random_state)
        if device == 'cuda':
            torch.cuda.manual_seed(random_state)
        
        logger.info("Running Variational Autoencoder anomaly detection")
        
        # Normalize data
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(embeddings)
        
        # Convert to tensor dataset
        tensor_data = torch.tensor(scaled_data, dtype=torch.float32)
        dataset = data_utils.TensorDataset(tensor_data)
        dataloader = data_utils.DataLoader(dataset, batch_size=batch_size, shuffle=True)
        
        # Set up VAE model
        input_dim = embeddings.shape[1]
        
        # Try to load pre-trained model if requested
        model = None
        if load_model_path and os.path.exists(load_model_path):
            try:
                logger.info(f"Loading pre-trained VAE model from {load_model_path}")
                model = torch.load(load_model_path, map_location=device)
            except Exception as e:
                logger.warning(f"Error loading pre-trained model: {e}")
        
        # Create new model if loading failed or wasn't requested
        if model is None:
            model = setup_vae_model(input_dim, latent_dim, device)
            
            # Train VAE
            model = train_vae_model(model, dataloader, epochs, lr, device)
        
        # Compute reconstruction error and KL divergence for anomaly scoring
        model.eval()
        anomaly_scores = []
        
        with torch.no_grad():
            for batch in dataloader:
                batch_x = batch[0].to(device)
                recon_x, mu, logvar = model(batch_x)
                
                # Compute reconstruction error
                recon_error = torch.mean(torch.square(batch_x - recon_x), dim=1)
                
                # Compute KL divergence
                kl_div = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp(), dim=1)
                
                # Combined score
                score = recon_error + 0.1 * kl_div  # Weight KL less to focus on reconstruction
                
                anomaly_scores.append(score.cpu())
        
        # Concatenate all scores
        all_scores = torch.cat(anomaly_scores).numpy()
        
        # Determine threshold based on contamination
        threshold = np.percentile(all_scores, 100 * (1 - contamination))
        labels = np.where(all_scores >= threshold, 1, 0)
        
        return labels, all_scores, model
        
    except ImportError:
        logger.warning("TensorFlow or PyTorch not found. Skipping VAE anomaly detection.")
        return np.zeros(len(embeddings)), np.zeros(len(embeddings)), None

def deep_svdd_detection(embeddings, contamination=0.05, random_state=42,
                      latent_dim=64, epochs=100, lr=0.001, batch_size=32, device='cpu',
                      load_model_path=None):
    """Detect anomalies using Deep SVDD."""
    try:
        import torch
        import torch.utils.data as data_utils
        
        # Set random seeds for reproducibility
        torch.manual_seed(random_state)
        if device == 'cuda':
            torch.cuda.manual_seed(random_state)
        
        logger.info("Running Deep SVDD anomaly detection")
        
        # Normalize data
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(embeddings)
        
        # Convert to tensor dataset
        tensor_data = torch.tensor(scaled_data, dtype=torch.float32)
        dataset = data_utils.TensorDataset(tensor_data)
        dataloader = data_utils.DataLoader(dataset, batch_size=batch_size, shuffle=True)
        
        # Set up Deep SVDD model
        input_dim = embeddings.shape[1]
        
        # Try to load pre-trained model if requested
        model = None
        if load_model_path and os.path.exists(load_model_path):
            try:
                logger.info(f"Loading pre-trained DeepSVDD model from {load_model_path}")
                model = torch.load(load_model_path, map_location=device)
            except Exception as e:
                logger.warning(f"Error loading pre-trained model: {e}")
        
        # Create new model if loading failed or wasn't requested
        if model is None:
            model = setup_deep_svdd_model(input_dim, latent_dim, device)
            
            # Train Deep SVDD
            model = train_deep_svdd_model(model, dataloader, epochs, lr, device)
        
        # Compute distances to center for anomaly scoring
        model.eval()
        anomaly_scores = []
        
        with torch.no_grad():
            for batch in dataloader:
                batch_x = batch[0].to(device)
                output = model(batch_x)
                
                # Compute distance to center
                dist = torch.sum((output - model.center) ** 2, dim=1)
                anomaly_scores.append(dist.cpu())
        
        # Concatenate all scores
        all_scores = torch.cat(anomaly_scores).numpy()
        
        # Determine threshold based on contamination
        threshold = np.percentile(all_scores, 100 * (1 - contamination))
        labels = np.where(all_scores >= threshold, 1, 0)
        
        return labels, all_scores, model
        
    except ImportError:
        logger.warning("PyTorch not found. Skipping Deep SVDD anomaly detection.")
        return np.zeros(len(embeddings)), np.zeros(len(embeddings)), None

def plot_anomalies(embeddings_2d, labels, scores, method_name, dim_red_method, output_dir, ids=None, annotations=None):
    """Plot anomalies in 2D embedding space."""
    plt.figure(figsize=(12, 10))
    
    # Create a colormap for the anomaly scores
    normal_mask = labels == 0
    anomaly_mask = labels == 1
    
    # Plot normal points
    plt.scatter(
        embeddings_2d[normal_mask, 0],
        embeddings_2d[normal_mask, 1],
        c='blue',
        alpha=0.5,
        label='Normal'
    )
    
    # Plot anomalies with color intensity based on score
    if np.any(anomaly_mask):
        anomaly_colors = plt.cm.OrRd(scores[anomaly_mask] / max(scores[anomaly_mask]))
        scatter = plt.scatter(
            embeddings_2d[anomaly_mask, 0],
            embeddings_2d[anomaly_mask, 1],
            c=scores[anomaly_mask],
            cmap='OrRd',
            label='Anomalies'
        )
        plt.colorbar(scatter, label='Anomaly Score')
        
        # Add labels for the top anomalies
        if ids is not None:
            top_n = min(10, sum(anomaly_mask))
            top_indices = np.argsort(scores)[-top_n:]
            for idx in top_indices:
                if labels[idx] == 1:  # Only label anomalies
                    plt.annotate(
                        ids[idx],
                        (embeddings_2d[idx, 0], embeddings_2d[idx, 1]),
                        fontsize=8,
                        alpha=0.7
                    )
    
    plt.title(f'Anomaly Detection using {method_name}\n2D Embedding with {dim_red_method.upper()}')
    plt.legend()
    plt.tight_layout()
    
    # Save the plot
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, f'{method_name}_{dim_red_method}.png'), dpi=300)
    plt.close()
    
    # If annotations are provided, create a plot colored by taxonomy
    if annotations is not None:
        plt.figure(figsize=(14, 12))
        
        # Get unique taxonomy groups
        unique_taxa = set(annotations.values())
        
        # Remove None or empty values
        unique_taxa = [t for t in unique_taxa if t and t != 'unknown']
        
        # Limit to top 10 most common taxa + unknown + anomaly
        if len(unique_taxa) > 10:
            taxa_counts = Counter([annotations.get(id, 'unknown') for id in ids])
            top_taxa = [t for t, _ in taxa_counts.most_common(10) if t and t != 'unknown']
            if 'unknown' not in top_taxa:
                top_taxa.append('unknown')
        else:
            top_taxa = unique_taxa
            if 'unknown' not in top_taxa:
                top_taxa.append('unknown')
                
        # Create a colormap for taxonomies
        colors = plt.cm.tab20(np.linspace(0, 1, len(top_taxa)))
        taxa_color_map = {taxon: color for taxon, color in zip(top_taxa, colors)}
        
        # Add "Other" category
        taxa_color_map['Other'] = plt.cm.tab20(0.95)
        
        # Plot points colored by taxonomy
        for i, id in enumerate(ids):
            taxon = annotations.get(id, 'unknown')
            
            if taxon not in top_taxa and taxon != 'unknown':
                taxon = 'Other'
                
            # Determine if this point is an anomaly
            is_anomaly = labels[i] == 1
            
            # Adjust marker properties based on anomaly status
            marker_size = 50 if is_anomaly else 20
            marker_edge_color = 'red' if is_anomaly else None
            marker_line_width = 1 if is_anomaly else 0
            
            plt.scatter(
                embeddings_2d[i, 0],
                embeddings_2d[i, 1],
                c=[taxa_color_map.get(taxon, taxa_color_map.get('unknown'))],
                s=marker_size,
                edgecolors=marker_edge_color,
                linewidths=marker_line_width,
                alpha=0.7
            )
        
        # Create legend patches for taxonomy
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=taxa_color_map[taxon], label=taxon)
            for taxon in taxa_color_map
        ]
        
        # Add legend for anomalies
        legend_elements.append(
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue',
                       markeredgecolor='red', markersize=10, label='Anomaly')
        )
        
        plt.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.15, 1))
        plt.title(f'Taxonomy and Anomalies ({method_name}, {dim_red_method.upper()})')
        plt.tight_layout()
        
        # Save the taxonomy plot
        plt.savefig(os.path.join(output_dir, f'{method_name}_{dim_red_method}_taxonomy.png'), dpi=300)
        plt.close()

def ensemble_anomaly_detection(embeddings, methods, contamination=0.05, ensemble_threshold=0.5, 
                              random_state=42, threads=8, latent_dim=64, epochs=100, lr=0.001, 
                              batch_size=32, device='cpu', save_dir=None, load_model=False):
    """Run multiple anomaly detection methods and combine results."""
    method_labels = {}
    method_scores = {}
    
    # Run each method
    for method in methods:
        if method == 'isolation_forest':
            labels, scores, _ = isolation_forest_detection(
                embeddings, contamination, random_state, threads
            )
        elif method == 'dbscan':
            labels, scores, _ = dbscan_detection(
                embeddings, contamination
            )
        elif method == 'lof':
            labels, scores, _ = local_outlier_factor_detection(
                embeddings, contamination, threads
            )
        elif method == 'one_class_svm':
            labels, scores, _ = one_class_svm_detection(
                embeddings, contamination
            )
        elif method == 'vae':
            load_model_path = None
            if load_model and save_dir:
                load_model_path = os.path.join(save_dir, f'vae_model.pt')
            
            labels, scores, model = variational_autoencoder_detection(
                embeddings, contamination, random_state, latent_dim, 
                epochs, lr, batch_size, device, load_model_path
            )
            
            # Save model if successfully trained
            if model is not None and save_dir:
                os.makedirs(save_dir, exist_ok=True)
                torch.save(model, os.path.join(save_dir, f'vae_model.pt'))
                
        elif method == 'deep_svdd':
            load_model_path = None
            if load_model and save_dir:
                load_model_path = os.path.join(save_dir, f'deep_svdd_model.pt')
            
            labels, scores, model = deep_svdd_detection(
                embeddings, contamination, random_state, latent_dim, 
                epochs, lr, batch_size, device, load_model_path
            )
            
            # Save model if successfully trained
            if model is not None and save_dir:
                os.makedirs(save_dir, exist_ok=True)
                torch.save(model, os.path.join(save_dir, f'deep_svdd_model.pt'))
        else:
            logger.warning(f"Unknown method: {method}, skipping")
            continue
        
        method_labels[method] = labels
        method_scores[method] = scores
    
    # Combine results
    if not method_labels:  # No methods ran successfully
        return np.zeros(len(embeddings)), np.zeros(len(embeddings))
    
    # Create a matrix of all labels
    all_labels = np.vstack([method_labels[m] for m in method_labels])
    
    # Vote (a point is an anomaly if at least ensemble_threshold of methods agree)
    vote_threshold = max(1, int(ensemble_threshold * len(method_labels)))
    ensemble_labels = np.sum(all_labels, axis=0) >= vote_threshold
    
    # Compute ensemble scores by averaging normalized scores
    all_scores = []
    for method in method_scores:
        # Normalize scores to [0, 1]
        method_score = method_scores[method]
        if np.max(method_score) > np.min(method_score):
            norm_score = (method_score - np.min(method_score)) / (np.max(method_score) - np.min(method_score))
        else:
            norm_score = np.zeros_like(method_score)
        all_scores.append(norm_score)
    
    ensemble_scores = np.mean(all_scores, axis=0)
    
    return ensemble_labels.astype(int), ensemble_scores

def incorporate_precomputed_anomalies(precomputed_anomalies, ids, ensemble_results=None):
    """Incorporate pre-computed anomalies from the embedding step."""
    if precomputed_anomalies is None:
        return ensemble_results
    
    # Extract labels and scores from precomputed anomalies
    anomaly_labels = precomputed_anomalies['labels']
    anomaly_scores = precomputed_anomalies['scores']
    
    # Create arrays aligned with the ids
    precomputed_labels = np.array([1 if anomaly_labels.get(id, False) else 0 for id in ids])
    precomputed_scores = np.array([anomaly_scores.get(id, 0.0) for id in ids])
    
    # If no ensemble results, return precomputed anomalies
    if ensemble_results is None:
        return precomputed_labels, precomputed_scores
    
    # Otherwise, combine with ensemble results
    ensemble_labels, ensemble_scores = ensemble_results
    
    # Use OR operation for labels (anomalous if either method says so)
    combined_labels = np.logical_or(precomputed_labels, ensemble_labels).astype(int)
    
    # For scores, take the maximum of the two
    # Normalize scores to [0, 1] first
    if np.max(ensemble_scores) > np.min(ensemble_scores):
        norm_ensemble = (ensemble_scores - np.min(ensemble_scores)) / (np.max(ensemble_scores) - np.min(ensemble_scores))
    else:
        norm_ensemble = np.zeros_like(ensemble_scores)
    
    if np.max(precomputed_scores) > np.min(precomputed_scores):
        norm_precomputed = (precomputed_scores - np.min(precomputed_scores)) / (np.max(precomputed_scores) - np.min(precomputed_scores))
    else:
        norm_precomputed = np.zeros_like(precomputed_scores)
    
    combined_scores = np.maximum(norm_ensemble, norm_precomputed)
    
    return combined_labels, combined_scores

def run_anomaly_detection(embeddings, ids, annotations, methods, save_dir, visualize_dir,
                         contamination=0.05, dim_reduction_methods=None, data_type='protein',
                         perplexity=30, random_state=42, n_jobs=8, 
                         precomputed_anomalies=None, ensemble=False, ensemble_threshold=0.5,
                         latent_dim=64, epochs=100, learning_rate=0.001, batch_size=32,
                         device='cpu', load_models=False):
    """Run multiple anomaly detection methods and dimensionality reduction techniques."""
    results = []
    
    # Apply dimensionality reduction
    if dim_reduction_methods is None:
        dim_reduction_methods = ['umap', 'tsne']
    
    dim_reduced, scaler = reduce_dimensions(
        embeddings, 
        dim_reduction_methods,
        perplexity=perplexity,
        random_state=random_state,
        n_jobs=n_jobs
    )
    
    # Convert annotations dict to array
    annotation_array = np.array([annotations.get(id, 'unknown') for id in ids]) if annotations else None
    
    # Run anomaly detection methods
    if ensemble:
        logger.info(f"Running ensemble anomaly detection with methods: {', '.join(methods)}")
        
        # Run ensemble
        ensemble_model_dir = os.path.join(save_dir, f"{data_type}_ensemble") if save_dir else None
        labels, scores = ensemble_anomaly_detection(
            embeddings, methods, contamination, ensemble_threshold,
            random_state, n_jobs, latent_dim, epochs, learning_rate,
            batch_size, device, ensemble_model_dir, load_models
        )
        
        # Incorporate precomputed anomalies if available
        if precomputed_anomalies:
            labels, scores = incorporate_precomputed_anomalies(
                precomputed_anomalies, ids, (labels, scores)
            )
        
        # Create visualizations
        if visualize_dir:
            os.makedirs(visualize_dir, exist_ok=True)
            
            for dim_method, embedding_2d in dim_reduced.items():
                plot_anomalies(
                    embedding_2d,
                    labels,
                    scores,
                    "Ensemble",
                    dim_method,
                    visualize_dir,
                    ids,
                    annotations
                )
        
        # Store ensemble results
        for i, (id, label, score) in enumerate(zip(ids, labels, scores)):
            results.append({
                'id': id,
                'method': 'ensemble',
                'anomaly': bool(label),
                'score': score,
                'taxonomy': annotations.get(id, 'unknown') if annotations else 'unknown'
            })
    
    else:
        # Run each anomaly detection method separately
        for method in methods:
            logger.info(f"Running {method} anomaly detection")
            
            if method == 'isolation_forest':
                labels, scores, model = isolation_forest_detection(
                    embeddings, contamination, random_state, n_jobs
                )
            elif method == 'dbscan':
                labels, scores, model = dbscan_detection(
                    embeddings, contamination
                )
            elif method == 'lof':
                labels, scores, model = local_outlier_factor_detection(
                    embeddings, contamination, n_jobs
                )
            elif method == 'one_class_svm':
                labels, scores, model = one_class_svm_detection(
                    embeddings, contamination
                )
            elif method == 'vae':
                # Set up model path for loading/saving
                load_model_path = None
                if load_models and save_dir:
                    load_model_path = os.path.join(save_dir, f"{data_type}_vae_model.pt")
                
                labels, scores, model = variational_autoencoder_detection(
                    embeddings, contamination, random_state, latent_dim, 
                    epochs, learning_rate, batch_size, device, load_model_path
                )
                
                # Save model if successfully trained
                if model is not None and save_dir:
                    os.makedirs(save_dir, exist_ok=True)
                    torch.save(model, os.path.join(save_dir, f"{data_type}_vae_model.pt"))
            
            elif method == 'deep_svdd':
                # Set up model path for loading/saving
                load_model_path = None
                if load_models and save_dir:
                    load_model_path = os.path.join(save_dir, f"{data_type}_deep_svdd_model.pt")
                
                labels, scores, model = deep_svdd_detection(
                    embeddings, contamination, random_state, latent_dim, 
                    epochs, learning_rate, batch_size, device, load_model_path
                )
                
                # Save model if successfully trained
                if model is not None and save_dir:
                    os.makedirs(save_dir, exist_ok=True)
                    torch.save(model, os.path.join(save_dir, f"{data_type}_deep_svdd_model.pt"))
            
            else:
                logger.warning(f"Unknown method: {method}, skipping")
                continue
            
            # Incorporate precomputed anomalies for this method if available
            if precomputed_anomalies and precomputed_anomalies['method'] == method:
                labels, scores = incorporate_precomputed_anomalies(
                    precomputed_anomalies, ids, (labels, scores)
                )
            
            # Save model if available
            if model is not None and save_dir and method not in ['vae', 'deep_svdd']:  # These are saved separately above
                os.makedirs(save_dir, exist_ok=True)
                try:
                    dump(model, os.path.join(save_dir, f'{data_type}_{method}_model.joblib'))
                    dump(scaler, os.path.join(save_dir, f'{data_type}_scaler.joblib'))
                except Exception as e:
                    logger.warning(f"Error saving model: {e}")
            
            # Create visualizations
            if visualize_dir:
                os.makedirs(visualize_dir, exist_ok=True)
                
                for dim_method, embedding_2d in dim_reduced.items():
                    plot_anomalies(
                        embedding_2d,
                        labels,
                        scores,
                        method,
                        dim_method,
                        visualize_dir,
                        ids,
                        annotations
                    )
            
            # Store results
            for i, (id, label, score) in enumerate(zip(ids, labels, scores)):
                results.append({
                    'id': id,
                    'method': method,
                    'anomaly': bool(label),
                    'score': score,
                    'taxonomy': annotations.get(id, 'unknown') if annotations else 'unknown'
                })
    
    # Convert results to DataFrame
    results_df = pd.DataFrame(results)
    
    return results_df

def main():
    args = parse_args()
    
    # Set number of threads
    os.environ["OMP_NUM_THREADS"] = str(args.threads)
    
    # Set up device for PyTorch
    device = 'cuda' if torch.cuda.is_available() and args.threads > 0 else 'cpu'
    
    # Convert method strings to lists
    detection_methods = [m.strip() for m in args.methods.split(',')]
    dim_reduction_methods = [m.strip() for m in args.dimensions_reduction.split(',')]
    
    # Process protein embeddings if provided
    if args.protein_embeddings and args.protein_ids:
        logger.info("Processing protein embeddings")
        
        # Load protein data
        protein_embeddings, protein_ids, protein_annotations, precomputed_protein_anomalies = load_data(
            args.protein_embeddings, args.protein_ids, args.protein_annotations, args.protein_anomalies
        )
        
        # Set up save directory for protein models
        protein_save_dir = os.path.join(args.save_models, "protein") if args.save_models else None
        
        # Run anomaly detection
        protein_results = run_anomaly_detection(
            protein_embeddings,
            protein_ids,
            protein_annotations,
            detection_methods,
            protein_save_dir,
            args.visualize,
            args.contamination,
            dim_reduction_methods,
            'protein',
            args.perplexity,
            args.random_state,
            args.threads,
            precomputed_protein_anomalies,
            args.ensemble,
            args.ensemble_threshold,
            args.latent_dim,
            args.epochs,
            args.learning_rate,
            min(32, len(protein_ids)),  # Adjust batch size to not exceed dataset size
            device,
            args.load_models
        )
        
        # Save results
        if args.output_protein:
            protein_results.to_csv(args.output_protein, sep='\t', index=False)
    
    # Process DNA embeddings if provided
    if args.dna_embeddings and args.contig_ids:
        logger.info("Processing DNA embeddings")
        
        # Load DNA data
        dna_embeddings, contig_ids, contig_annotations, precomputed_dna_anomalies = load_data(
            args.dna_embeddings, args.contig_ids, args.contig_annotations, args.dna_anomalies
        )
        
        # Set up save directory for DNA models
        dna_save_dir = os.path.join(args.save_models, "dna") if args.save_models else None
        
        # Run anomaly detection
        dna_results = run_anomaly_detection(
            dna_embeddings,
            contig_ids,
            contig_annotations,
            detection_methods,
            dna_save_dir,
            args.visualize,
            args.contamination,
            dim_reduction_methods,
            'dna',
            args.perplexity,
            args.random_state,
            args.threads,
            precomputed_dna_anomalies,
            args.ensemble,
            args.ensemble_threshold,
            args.latent_dim,
            args.epochs,
            args.learning_rate,
            min(32, len(contig_ids)),  # Adjust batch size to not exceed dataset size
            device,
            args.load_models
        )
        
        # Save results
        if args.output_dna:
            dna_results.to_csv(args.output_dna, sep='\t', index=False)
    
    logger.info("Anomaly detection completed successfully")

if __name__ == "__main__":
    main()
