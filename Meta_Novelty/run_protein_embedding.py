#!/usr/bin/env python
"""
Generate protein embeddings from FASTA sequences using various embedding models.
Models supported: ESM2, ProtT5, ESM1b
Anomaly detection methods: Isolation Forest, VAE, DeepSVDD
"""

import argparse
import numpy as np
import os
import sys
import torch
from Bio import SeqIO
from tqdm import tqdm
import logging
import warnings
import json
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import IsolationForest
warnings.filterwarnings("ignore")

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('protein_embedding')

def parse_args():
    parser = argparse.ArgumentParser(description='Generate protein embeddings')
    parser.add_argument('--input', required=True, help='Input protein FASTA file')
    parser.add_argument('--output', required=True, help='Output NPY file for embeddings')
    parser.add_argument('--ids', required=True, help='Output text file for protein IDs')
    parser.add_argument('--model', default='esm2', choices=['esm2', 'prot_t5', 'esm1b', 'vae'], 
                        help='Embedding model to use')
    parser.add_argument('--batch_size', type=int, default=32, 
                        help='Batch size for embedding generation')
    parser.add_argument('--max_seq_len', type=int, default=1024, 
                        help='Maximum sequence length')
    parser.add_argument('--threads', type=int, default=8, 
                        help='Number of CPU threads')
    parser.add_argument('--device', default='cpu', choices=['cpu', 'cuda'],
                        help='Device to use for computation')
    
    # Anomaly detection parameters
    parser.add_argument('--detect_anomalies', action='store_true',
                        help='Enable anomaly detection during embedding generation')
    parser.add_argument('--anomaly_method', default='isolation_forest', 
                        choices=['isolation_forest', 'vae', 'deep_svdd'],
                        help='Anomaly detection method to use')
    parser.add_argument('--contamination', type=float, default=0.05,
                        help='Expected fraction of anomalies in the data')
    parser.add_argument('--filter_anomalies', action='store_true',
                        help='Filter out anomalies from the embedding output')
    parser.add_argument('--save_anomalies', action='store_true',
                        help='Save anomalies to a separate file')
    parser.add_argument('--anomaly_output', help='Output file for anomaly results (JSON)')
    
    # VAE/DeepSVDD specific parameters
    parser.add_argument('--latent_dim', type=int, default=64,
                        help='Latent dimension size for VAE or DeepSVDD')
    parser.add_argument('--epochs', type=int, default=100,
                        help='Number of training epochs for deep learning models')
    parser.add_argument('--learning_rate', type=float, default=0.001,
                        help='Learning rate for neural network training')
    
    return parser.parse_args()

def load_sequences(fasta_file):
    """Load protein sequences from FASTA file."""
    logger.info(f"Loading sequences from {fasta_file}")
    sequences = []
    ids = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
        ids.append(record.id)
    return sequences, ids

def setup_esm2_model(device, max_seq_len):
    """Set up ESM2 model."""
    try:
        import esm
        logger.info("Loading ESM2 model")
        # Load ESM2 model (35M parameter version for efficiency)
        model, alphabet = esm.pretrained.esm2_t12_35M_UR50D()
        model = model.to(device)
        model.eval()
        batch_converter = alphabet.get_batch_converter(max_seq_len=max_seq_len)
        return model, batch_converter
    except ImportError:
        logger.error("ESM module not found. Install with: pip install fair-esm")
        sys.exit(1)

def setup_prot_t5_model(device, max_seq_len):
    """Set up ProtT5 model."""
    try:
        from transformers import T5EncoderModel, T5Tokenizer
        logger.info("Loading ProtT5 model")
        tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_uniref50", do_lower_case=False)
        model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50")
        model = model.to(device)
        model.eval()
        return model, tokenizer
    except ImportError:
        logger.error("Transformers module not found. Install with: pip install transformers")
        sys.exit(1)

def setup_esm1b_model(device, max_seq_len):
    """Set up ESM1b model."""
    try:
        import esm
        logger.info("Loading ESM1b model")
        model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
        model = model.to(device)
        model.eval()
        batch_converter = alphabet.get_batch_converter(max_seq_len=max_seq_len)
        return model, batch_converter
    except ImportError:
        logger.error("ESM module not found. Install with: pip install fair-esm")
        sys.exit(1)

def setup_vae_model(input_dim, latent_dim, device):
    """Set up VAE model for embedding generation and anomaly detection."""
    try:
        import torch.nn as nn
        import torch.nn.functional as F
        
        class ProteinVAE(nn.Module):
            def __init__(self, input_dim, hidden_dim=512, latent_dim=64):
                super(ProteinVAE, self).__init__()
                
                # Encoder
                self.fc1 = nn.Linear(input_dim, hidden_dim)
                self.fc2 = nn.Linear(hidden_dim, hidden_dim)
                self.fc31 = nn.Linear(hidden_dim, latent_dim)  # mean
                self.fc32 = nn.Linear(hidden_dim, latent_dim)  # log variance
                
                # Decoder
                self.fc4 = nn.Linear(latent_dim, hidden_dim)
                self.fc5 = nn.Linear(hidden_dim, hidden_dim)
                self.fc6 = nn.Linear(hidden_dim, input_dim)
                
                self.latent_dim = latent_dim
            
            def encode(self, x):
                h1 = F.relu(self.fc1(x))
                h2 = F.relu(self.fc2(h1))
                return self.fc31(h2), self.fc32(h2)
            
            def reparameterize(self, mu, logvar):
                std = torch.exp(0.5 * logvar)
                eps = torch.randn_like(std)
                return mu + eps * std
            
            def decode(self, z):
                h4 = F.relu(self.fc4(z))
                h5 = F.relu(self.fc5(h4))
                return self.fc6(h5)
            
            def forward(self, x):
                mu, logvar = self.encode(x)
                z = self.reparameterize(mu, logvar)
                return self.decode(z), mu, logvar
            
            def get_embedding(self, x):
                mu, _ = self.encode(x)
                return mu
            
            def compute_anomaly_score(self, x):
                x_recon, mu, logvar = self.forward(x)
                # Reconstruction error
                recon_error = F.mse_loss(x_recon, x, reduction='none').sum(dim=1)
                # KL divergence
                kl_div = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp(), dim=1)
                # Combined anomaly score (weighted sum)
                anomaly_score = recon_error + kl_div
                return anomaly_score
        
        logger.info(f"Creating VAE model with input dim {input_dim} and latent dim {latent_dim}")
        model = ProteinVAE(input_dim, hidden_dim=512, latent_dim=latent_dim)
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
            def __init__(self, input_dim, hidden_dim=512, output_dim=64):
                super(DeepSVDD, self).__init__()
                self.fc1 = nn.Linear(input_dim, hidden_dim)
                self.bn1 = nn.BatchNorm1d(hidden_dim)
                self.fc2 = nn.Linear(hidden_dim, hidden_dim)
                self.bn2 = nn.BatchNorm1d(hidden_dim)
                self.fc3 = nn.Linear(hidden_dim, output_dim)
                
                # Initialize center to None (will be set during training)
                self.center = None
                self.radius = 0.0
                
            def forward(self, x):
                x = F.relu(self.bn1(self.fc1(x)))
                x = F.relu(self.bn2(self.fc2(x)))
                x = self.fc3(x)
                return x
            
            def get_embedding(self, x):
                return self.forward(x)
            
            def init_center(self, embeddings, eps=0.1):
                """Initialize hypersphere center as the mean of all embeddings."""
                self.center = torch.mean(embeddings, dim=0)
                # Add small constant to avoid collapse to a single point
                self.center[(abs(self.center) < eps) & (self.center < 0)] = -eps
                self.center[(abs(self.center) < eps) & (self.center >= 0)] = eps
            
            def compute_anomaly_score(self, x):
                embeddings = self.forward(x)
                if self.center is None:
                    raise ValueError("Center not initialized. Call init_center first.")
                    
                # Distance to center
                dist = torch.sum((embeddings - self.center) ** 2, dim=1)
                return dist
        
        logger.info(f"Creating DeepSVDD model with input dim {input_dim} and output dim {latent_dim}")
        model = DeepSVDD(input_dim, hidden_dim=512, output_dim=latent_dim)
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
    for epoch in tqdm(range(epochs), desc="Training VAE"):
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
        
        # Log progress
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
    for epoch in tqdm(range(epochs), desc="Training DeepSVDD"):
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
        
        # Log progress
        if (epoch + 1) % 10 == 0:
            logger.info(f"Epoch {epoch+1}/{epochs}, Loss: {total_loss / len(dataloader.dataset):.4f}")
    
    # Compute radius (95th percentile of distances)
    model.eval()
    distances = []
    with torch.no_grad():
        for batch_x in dataloader:
            batch_x = batch_x.to(device)
            embeddings = model(batch_x)
            dist = torch.sum((embeddings - model.center) ** 2, dim=1)
            distances.append(dist)
    
    all_distances = torch.cat(distances, dim=0)
    model.radius = torch.quantile(all_distances, 0.95).item()
    logger.info(f"DeepSVDD radius set to {model.radius:.4f}")
    
    return model

def compute_kmers(sequences, k=3):
    """Compute k-mer frequencies for protein sequences."""
    from collections import Counter
    import itertools
    
    # Generate all possible k-mers for amino acids
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    all_kmers = [''.join(p) for p in itertools.product(amino_acids, repeat=k)]
    kmer_to_idx = {kmer: i for i, kmer in enumerate(all_kmers)}
    
    # Compute k-mer frequencies for each sequence
    kmer_features = []
    
    for seq in tqdm(sequences, desc=f"Computing {k}-mer frequencies"):
        kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
        counter = Counter(kmers)
        
        # Create feature vector
        feature_vec = np.zeros(len(kmer_to_idx))
        for kmer, count in counter.items():
            if kmer in kmer_to_idx:
                feature_vec[kmer_to_idx[kmer]] = count / len(kmers)
        
        kmer_features.append(feature_vec)
    
    return np.array(kmer_features)

def generate_esm2_embeddings(model, batch_converter, sequences, ids, batch_size, device):
    """Generate embeddings using ESM2 model."""
    embeddings = {}
    
    for i in tqdm(range(0, len(sequences), batch_size)):
        batch_ids = ids[i:i+batch_size]
        batch_seqs = sequences[i:i+batch_size]
        
        batch_data = [(id, seq) for id, seq in zip(batch_ids, batch_seqs)]
        
        with torch.no_grad():
            batch_labels, batch_strs, batch_tokens = batch_converter(batch_data)
            batch_tokens = batch_tokens.to(device)
            results = model(batch_tokens, repr_layers=[12])
            
            # Extract per-protein embeddings (average over sequence length, excluding special tokens)
            token_embeddings = results["representations"][12]
            
            for i, id in enumerate(batch_ids):
                # Exclude special tokens (first and last tokens)
                seq_len = len(batch_seqs[i])
                emb = token_embeddings[i, 1:seq_len+1].mean(0).cpu().numpy()
                embeddings[id] = emb
    
    return embeddings

def generate_prot_t5_embeddings(model, tokenizer, sequences, ids, batch_size, device, max_seq_len):
    """Generate embeddings using ProtT5 model."""
    embeddings = {}
    
    for i in tqdm(range(0, len(sequences), batch_size)):
        batch_ids = ids[i:i+batch_size]
        batch_seqs = sequences[i:i+batch_size]
        
        # Reformat sequences for ProtT5 (space between each amino acid)
        batch_seqs = [" ".join(list(seq[:max_seq_len])) for seq in batch_seqs]
        
        with torch.no_grad():
            # Tokenize and generate embeddings
            encoded_input = tokenizer(batch_seqs, return_tensors="pt", padding=True, truncation=True)
            encoded_input = {k: v.to(device) for k, v in encoded_input.items()}
            outputs = model(**encoded_input)
            
            # Mean over sequence length (exclude padding)
            for j, id in enumerate(batch_ids):
                # Get attention mask to exclude padding tokens
                attention_mask = encoded_input["attention_mask"][j]
                embeddings_sum = torch.sum(outputs.last_hidden_state[j] * attention_mask.unsqueeze(-1), dim=0)
                embeddings_avg = embeddings_sum / attention_mask.sum()
                embeddings[id] = embeddings_avg.cpu().numpy()
    
    return embeddings

def generate_esm1b_embeddings(model, batch_converter, sequences, ids, batch_size, device):
    """Generate embeddings using ESM1b model."""
    embeddings = {}
    
    for i in tqdm(range(0, len(sequences), batch_size)):
        batch_ids = ids[i:i+batch_size]
        batch_seqs = sequences[i:i+batch_size]
        
        batch_data = [(id, seq) for id, seq in zip(batch_ids, batch_seqs)]
        
        with torch.no_grad():
            batch_labels, batch_strs, batch_tokens = batch_converter(batch_data)
            batch_tokens = batch_tokens.to(device)
            results = model(batch_tokens, repr_layers=[33])
            
            # Extract per-protein embeddings (average over sequence length, excluding special tokens)
            token_embeddings = results["representations"][33]
            
            for i, id in enumerate(batch_ids):
                # Exclude special tokens (first and last tokens)
                seq_len = len(batch_seqs[i])
                emb = token_embeddings[i, 1:seq_len+1].mean(0).cpu().numpy()
                embeddings[id] = emb
    
    return embeddings

def generate_vae_embeddings(input_features, ids, latent_dim, epochs, learning_rate, batch_size, device):
    """Generate embeddings using VAE model."""
    import torch
    import torch.utils.data as data_utils
    
    # Convert input features to torch dataset
    input_features_tensor = torch.tensor(input_features, dtype=torch.float32)
    dataset = data_utils.TensorDataset(input_features_tensor)
    dataloader = data_utils.DataLoader(dataset, batch_size=batch_size, shuffle=True)
    
    # Create and train VAE model
    input_dim = input_features.shape[1]
    model = setup_vae_model(input_dim, latent_dim, device)
    model = train_vae_model(model, dataloader, epochs, learning_rate, device)
    
    # Generate embeddings
    embeddings = {}
    anomaly_scores = {}
    
    model.eval()
    with torch.no_grad():
        for i in tqdm(range(0, len(input_features), batch_size), desc="Generating VAE embeddings"):
            batch_ids = ids[i:i+batch_size]
            batch_x = torch.tensor(input_features[i:i+batch_size], dtype=torch.float32).to(device)
            
            # Get embeddings and anomaly scores
            mu = model.get_embedding(batch_x)
            scores = model.compute_anomaly_score(batch_x)
            
            for j, id in enumerate(batch_ids):
                embeddings[id] = mu[j].cpu().numpy()
                anomaly_scores[id] = scores[j].item()
    
    return embeddings, anomaly_scores, model

def generate_deep_svdd_embeddings(input_features, ids, latent_dim, epochs, learning_rate, batch_size, device):
    """Generate embeddings using Deep SVDD model."""
    import torch
    import torch.utils.data as data_utils
    
    # Convert input features to torch dataset
    input_features_tensor = torch.tensor(input_features, dtype=torch.float32)
    dataset = data_utils.TensorDataset(input_features_tensor)
    dataloader = data_utils.DataLoader(dataset, batch_size=batch_size, shuffle=True)
    
    # Create and train Deep SVDD model
    input_dim = input_features.shape[1]
    model = setup_deep_svdd_model(input_dim, latent_dim, device)
    model = train_deep_svdd_model(model, dataloader, epochs, learning_rate, device)
    
    # Generate embeddings
    embeddings = {}
    anomaly_scores = {}
    
    model.eval()
    with torch.no_grad():
        for i in tqdm(range(0, len(input_features), batch_size), desc="Generating DeepSVDD embeddings"):
            batch_ids = ids[i:i+batch_size]
            batch_x = torch.tensor(input_features[i:i+batch_size], dtype=torch.float32).to(device)
            
            # Get embeddings and anomaly scores
            emb = model.get_embedding(batch_x)
            scores = model.compute_anomaly_score(batch_x)
            
            for j, id in enumerate(batch_ids):
                embeddings[id] = emb[j].cpu().numpy()
                anomaly_scores[id] = scores[j].item()
    
    return embeddings, anomaly_scores, model

def detect_anomalies_isolation_forest(embeddings_dict, contamination=0.05):
    """Detect anomalies using Isolation Forest."""
    # Extract embeddings and IDs
    ids = list(embeddings_dict.keys())
    embeddings_array = np.array([embeddings_dict[id] for id in ids])
    
    # Standardize data
    scaler = StandardScaler()
    embeddings_scaled = scaler.fit_transform(embeddings_array)
    
    # Apply Isolation Forest
    logger.info("Running Isolation Forest for anomaly detection")
    iso_forest = IsolationForest(contamination=contamination, random_state=42, n_jobs=-1)
    
    # Fit and predict
    y_pred = iso_forest.fit_predict(embeddings_scaled)
    
    # Convert to binary labels and scores
    # In Isolation Forest: -1 for anomalies, 1 for normal
    # We convert to: 1 for anomalies, 0 for normal
    anomaly_labels = np.where(y_pred == -1, 1, 0)
    anomaly_scores = -iso_forest.decision_function(embeddings_scaled)  # Higher score means more anomalous
    
    # Create result dictionaries
    anomaly_dict = {id: bool(label) for id, label in zip(ids, anomaly_labels)}
    score_dict = {id: float(score) for id, score in zip(ids, anomaly_scores)}
    
    return anomaly_dict, score_dict

def main():
    args = parse_args()
    
    # Set number of threads
    torch.set_num_threads(args.threads)
    os.environ["OMP_NUM_THREADS"] = str(args.threads)
    
    # Check if CUDA is available when requested
    if args.device == 'cuda' and not torch.cuda.is_available():
        logger.warning("CUDA requested but not available. Using CPU instead.")
        args.device = 'cpu'
    
    # Load sequences
    sequences, ids = load_sequences(args.input)
    logger.info(f"Loaded {len(sequences)} protein sequences")
    
    # Filter out sequences that are too long or too short
    filtered_seqs = []
    filtered_ids = []
    for seq, id in zip(sequences, ids):
        if len(seq) > 0 and len(seq) <= args.max_seq_len:
            filtered_seqs.append(seq)
            filtered_ids.append(id)
        else:
            logger.warning(f"Skipping sequence {id} with length {len(seq)} (too long or empty)")
    
    logger.info(f"After filtering: {len(filtered_seqs)} protein sequences")
    
    # Check if enough sequences for anomaly detection
    if args.detect_anomalies and len(filtered_seqs) < 10:
        logger.warning("Too few sequences for anomaly detection. Disabling.")
        args.detect_anomalies = False
    
    # For VAE model, we need to first compute k-mer features
    if args.model == 'vae':
        logger.info("Computing k-mer features for VAE model")
        kmer_features = compute_kmers(filtered_seqs, k=3)
        
        # Generate embeddings using VAE
        embeddings_dict, anomaly_scores, model = generate_vae_embeddings(
            kmer_features, filtered_ids, args.latent_dim, args.epochs, 
            args.learning_rate, args.batch_size, args.device
        )
        
        # Set anomaly labels based on scores and contamination
        threshold = np.percentile(list(anomaly_scores.values()), 100 * (1 - args.contamination))
        anomaly_labels = {id: score >= threshold for id, score in anomaly_scores.items()}
        
    else:
        # Set up the model based on user choice
        if args.model == 'esm2':
            model, batch_converter = setup_esm2_model(args.device, args.max_seq_len)
            embeddings_dict = generate_esm2_embeddings(
                model, batch_converter, filtered_seqs, filtered_ids, args.batch_size, args.device
            )
        elif args.model == 'prot_t5':
            model, tokenizer = setup_prot_t5_model(args.device, args.max_seq_len)
            embeddings_dict = generate_prot_t5_embeddings(
                model, tokenizer, filtered_seqs, filtered_ids, args.batch_size, args.device, args.max_seq_len
            )
        elif args.model == 'esm1b':
            model, batch_converter = setup_esm1b_model(args.device, args.max_seq_len)
            embeddings_dict = generate_esm1b_embeddings(
                model, batch_converter, filtered_seqs, filtered_ids, args.batch_size, args.device
            )
        
        # Detect anomalies if requested
        if args.detect_anomalies:
            if args.anomaly_method == 'isolation_forest':
                anomaly_labels, anomaly_scores = detect_anomalies_isolation_forest(
                    embeddings_dict, args.contamination
                )
            elif args.anomaly_method == 'deep_svdd':
                # Extract embeddings and prepare for DeepSVDD
                ids_list = list(embeddings_dict.keys())
                embeddings_array = np.array([embeddings_dict[id] for id in ids_list])
                
                # Apply DeepSVDD
                _, anomaly_scores, _ = generate_deep_svdd_embeddings(
                    embeddings_array, ids_list, args.latent_dim, args.epochs, 
                    args.learning_rate, args.batch_size, args.device
                )
                
                # Set anomaly labels based on scores and contamination
                threshold = np.percentile(list(anomaly_scores.values()), 100 * (1 - args.contamination))
                anomaly_labels = {id: score >= threshold for id, score in anomaly_scores.items()}
    
    # Filter embeddings if requested
    if args.detect_anomalies and args.filter_anomalies:
        logger.info("Filtering anomalies from embeddings")
        filtered_embeddings = {}
        filtered_ids_list = []
        
        for id in embeddings_dict:
            if not anomaly_labels[id]:  # Only keep normal samples
                filtered_embeddings[id] = embeddings_dict[id]
                filtered_ids_list.append(id)
        
        logger.info(f"After filtering anomalies: {len(filtered_embeddings)} embeddings")
        embeddings_dict = filtered_embeddings
        filtered_ids = filtered_ids_list
    
    # Save anomalies if requested
    if args.detect_anomalies and args.save_anomalies and args.anomaly_output:
        logger.info(f"Saving anomaly detection results to {args.anomaly_output}")
        anomaly_results = {
            'method': args.anomaly_method,
            'contamination': args.contamination,
            'num_samples': len(filtered_ids),
            'num_anomalies': sum(anomaly_labels.values()),
            'results': []
        }
        
        for id in filtered_ids:
            anomaly_results['results'].append({
                'id': id,
                'is_anomaly': anomaly_labels[id],
                'score': anomaly_scores[id]
            })
        
        # Save to JSON file
        os.makedirs(os.path.dirname(args.anomaly_output) or '.', exist_ok=True)
        with open(args.anomaly_output, 'w') as f:
            json.dump(anomaly_results, f, indent=2)
    
    # Convert dictionary to ordered arrays
    ordered_ids = list(embeddings_dict.keys())
    ordered_embeddings = np.array([embeddings_dict[id] for id in ordered_ids])
    
    # Save embeddings and IDs
    logger.info(f"Saving embeddings to {args.output}")
    np.save(args.output, ordered_embeddings)
    
    logger.info(f"Saving protein IDs to {args.ids}")
    with open(args.ids, 'w') as f:
        for id in ordered_ids:
            f.write(f"{id}\n")
    
    logger.info(f"Successfully generated embeddings for {len(ordered_ids)} proteins")

if __name__ == "__main__":
    main()#!/usr/bin/env python
"""
Generate protein embeddings from FASTA sequences using various embedding models.
Models supported: ESM2, ProtT5, ESM1b
"""

import argparse
import numpy as np
import os
import sys
import torch
from Bio import SeqIO
from tqdm import tqdm
import logging
import warnings
warnings.filterwarnings("ignore")

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('protein_embedding')

def parse_args():
    parser = argparse.ArgumentParser(description='Generate protein embeddings')
    parser.add_argument('--input', required=True, help='Input protein FASTA file')
    parser.add_argument('--output', required=True, help='Output NPY file for embeddings')
    parser.add_argument('--ids', required=True, help='Output text file for protein IDs')
    parser.add_argument('--model', default='esm2', choices=['esm2', 'prot_t5', 'esm1b'], 
                        help='Embedding model to use')
    parser.add_argument('--batch_size', type=int, default=32, 
                        help='Batch size for embedding generation')
    parser.add_argument('--max_seq_len', type=int, default=1024, 
                        help='Maximum sequence length')
    parser.add_argument('--threads', type=int, default=8, 
                        help='Number of CPU threads')
    parser.add_argument('--device', default='cpu', choices=['cpu', 'cuda'],
                        help='Device to use for computation')
    return parser.parse_args()

def load_sequences(fasta_file):
    """Load protein sequences from FASTA file."""
    logger.info(f"Loading sequences from {fasta_file}")
    sequences = []
    ids = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
        ids.append(record.id)
    return sequences, ids

def setup_esm2_model(device, max_seq_len):
    """Set up ESM2 model."""
    try:
        import esm
        logger.info("Loading ESM2 model")
        # Load ESM2 model (35M parameter version for efficiency)
        model, alphabet = esm.pretrained.esm2_t12_35M_UR50D()
        model = model.to(device)
        model.eval()
        batch_converter = alphabet.get_batch_converter(max_seq_len=max_seq_len)
        return model, batch_converter
    except ImportError:
        logger.error("ESM module not found. Install with: pip install fair-esm")
        sys.exit(1)

def setup_prot_t5_model(device, max_seq_len):
    """Set up ProtT5 model."""
    try:
        from transformers import T5EncoderModel, T5Tokenizer
        logger.info("Loading ProtT5 model")
        tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_uniref50", do_lower_case=False)
        model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50")
        model = model.to(device)
        model.eval()
        return model, tokenizer
    except ImportError:
        logger.error("Transformers module not found. Install with: pip install transformers")
        sys.exit(1)

def setup_esm1b_model(device, max_seq_len):
    """Set up ESM1b model."""
    try:
        import esm
        logger.info("Loading ESM1b model")
        model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
        model = model.to(device)
        model.eval()
        batch_converter = alphabet.get_batch_converter(max_seq_len=max_seq_len)
        return model, batch_converter
    except ImportError:
        logger.error("ESM module not found. Install with: pip install fair-esm")
        sys.exit(1)

def generate_esm2_embeddings(model, batch_converter, sequences, ids, batch_size, device):
    """Generate embeddings using ESM2 model."""
    embeddings = {}
    
    for i in tqdm(range(0, len(sequences), batch_size)):
        batch_ids = ids[i:i+batch_size]
        batch_seqs = sequences[i:i+batch_size]
        
        batch_data = [(id, seq) for id, seq in zip(batch_ids, batch_seqs)]
        
        with torch.no_grad():
            batch_labels, batch_strs, batch_tokens = batch_converter(batch_data)
            batch_tokens = batch_tokens.to(device)
            results = model(batch_tokens, repr_layers=[12])
            
            # Extract per-protein embeddings (average over sequence length, excluding special tokens)
            token_embeddings = results["representations"][12]
            
            for i, id in enumerate(batch_ids):
                # Exclude special tokens (first and last tokens)
                seq_len = len(batch_seqs[i])
                emb = token_embeddings[i, 1:seq_len+1].mean(0).cpu().numpy()
                embeddings[id] = emb
    
    return embeddings

def generate_prot_t5_embeddings(model, tokenizer, sequences, ids, batch_size, device, max_seq_len):
    """Generate embeddings using ProtT5 model."""
    embeddings = {}
    
    for i in tqdm(range(0, len(sequences), batch_size)):
        batch_ids = ids[i:i+batch_size]
        batch_seqs = sequences[i:i+batch_size]
        
        # Reformat sequences for ProtT5 (space between each amino acid)
        batch_seqs = [" ".join(list(seq[:max_seq_len])) for seq in batch_seqs]
        
        with torch.no_grad():
            # Tokenize and generate embeddings
            encoded_input = tokenizer(batch_seqs, return_tensors="pt", padding=True, truncation=True)
            encoded_input = {k: v.to(device) for k, v in encoded_input.items()}
            outputs = model(**encoded_input)
            
            # Mean over sequence length (exclude padding)
            for j, id in enumerate(batch_ids):
                # Get attention mask to exclude padding tokens
                attention_mask = encoded_input["attention_mask"][j]
                embeddings_sum = torch.sum(outputs.last_hidden_state[j] * attention_mask.unsqueeze(-1), dim=0)
                embeddings_avg = embeddings_sum / attention_mask.sum()
                embeddings[id] = embeddings_avg.cpu().numpy()
    
    return embeddings

def generate_esm1b_embeddings(model, batch_converter, sequences, ids, batch_size, device):
    """Generate embeddings using ESM1b model."""
    embeddings = {}
    
    for i in tqdm(range(0, len(sequences), batch_size)):
        batch_ids = ids[i:i+batch_size]
        batch_seqs = sequences[i:i+batch_size]
        
        batch_data = [(id, seq) for id, seq in zip(batch_ids, batch_seqs)]
        
        with torch.no_grad():
            batch_labels, batch_strs, batch_tokens = batch_converter(batch_data)
            batch_tokens = batch_tokens.to(device)
            results = model(batch_tokens, repr_layers=[33])
            
            # Extract per-protein embeddings (average over sequence length, excluding special tokens)
            token_embeddings = results["representations"][33]
            
            for i, id in enumerate(batch_ids):
                # Exclude special tokens (first and last tokens)
                seq_len = len(batch_seqs[i])
                emb = token_embeddings[i, 1:seq_len+1].mean(0).cpu().numpy()
                embeddings[id] = emb
    
    return embeddings

def main():
    args = parse_args()
    
    # Set number of threads
    torch.set_num_threads(args.threads)
    os.environ["OMP_NUM_THREADS"] = str(args.threads)
    
    # Check if CUDA is available when requested
    if args.device == 'cuda' and not torch.cuda.is_available():
        logger.warning("CUDA requested but not available. Using CPU instead.")
        args.device = 'cpu'
    
    # Load sequences
    sequences, ids = load_sequences(args.input)
    logger.info(f"Loaded {len(sequences)} protein sequences")
    
    # Filter out sequences that are too long or too short
    filtered_seqs = []
    filtered_ids = []
    for seq, id in zip(sequences, ids):
        if len(seq) > 0 and len(seq) <= args.max_seq_len:
            filtered_seqs.append(seq)
            filtered_ids.append(id)
        else:
            logger.warning(f"Skipping sequence {id} with length {len(seq)} (too long or empty)")
    
    logger.info(f"After filtering: {len(filtered_seqs)} protein sequences")
    
    # Set up the model based on user choice
    if args.model == 'esm2':
        model, batch_converter = setup_esm2_model(args.device, args.max_seq_len)
        embeddings_dict = generate_esm2_embeddings(
            model, batch_converter, filtered_seqs, filtered_ids, args.batch_size, args.device
        )
    elif args.model == 'prot_t5':
        model, tokenizer = setup_prot_t5_model(args.device, args.max_seq_len)
        embeddings_dict = generate_prot_t5_embeddings(
            model, tokenizer, filtered_seqs, filtered_ids, args.batch_size, args.device, args.max_seq_len
        )
    elif args.model == 'esm1b':
        model, batch_converter = setup_esm1b_model(args.device, args.max_seq_len)
        embeddings_dict = generate_esm1b_embeddings(
            model, batch_converter, filtered_seqs, filtered_ids, args.batch_size, args.device
        )
    
    # Convert dictionary to ordered arrays
    ordered_ids = list(embeddings_dict.keys())
    ordered_embeddings = np.array([embeddings_dict[id] for id in ordered_ids])
    
    # Save embeddings and IDs
    logger.info(f"Saving embeddings to {args.output}")
    np.save(args.output, ordered_embeddings)
    
    logger.info(f"Saving protein IDs to {args.ids}")
    with open(args.ids, 'w') as f:
        for id in ordered_ids:
            f.write(f"{id}\n")
    
    logger.info(f"Successfully generated embeddings for {len(ordered_ids)} proteins")

if __name__ == "__main__":
    main()

