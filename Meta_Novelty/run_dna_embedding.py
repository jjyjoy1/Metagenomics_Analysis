#!/usr/bin/env python
"""
Generate DNA embeddings from FASTA sequences using various embedding models.
Models supported: DNABERT, Nucleotide Transformer, K-mer based
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
from collections import defaultdict
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import IsolationForest
warnings.filterwarnings("ignore")

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('dna_embedding')

def parse_args():
    parser = argparse.ArgumentParser(description='Generate DNA embeddings')
    parser.add_argument('--input', required=True, help='Input DNA sequences file (FASTA or text)')
    parser.add_argument('--input_type', default='text', choices=['fasta', 'text'],
                        help='Input file type (FASTA or line-separated sequences)')
    parser.add_argument('--output', required=True, help='Output NPY file for embeddings')
    parser.add_argument('--ids', required=True, help='Output text file for sequence IDs')
    parser.add_argument('--model', default='dnabert', choices=['dnabert', 'nucleotide_transformer', 'k_mer', 'vae'],
                        help='Embedding model to use')
    parser.add_argument('--batch_size', type=int, default=8, 
                        help='Batch size for embedding generation')
    parser.add_argument('--max_seq_len', type=int, default=2048, 
                        help='Maximum sequence length')
    parser.add_argument('--kmer_size', type=int, default=6,
                        help='K-mer size for k-mer based embeddings')
    parser.add_argument('--threads', type=int, default=8, 
                        help='Number of CPU threads')
    parser.add_argument('--device', default='cpu', choices=['cpu', 'cuda'],
                        help='Device to use for computation')
    parser.add_argument('--sliding_window', type=int, default=500,
                        help='Sliding window size for long sequences')
    parser.add_argument('--sliding_step', type=int, default=250,
                        help='Step size for sliding window')
    
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

def load_sequences(file_path, input_type):
    """Load DNA sequences from file."""
    logger.info(f"Loading sequences from {file_path}")
    sequences = []
    ids = []
    
    if input_type == 'fasta':
        for record in SeqIO.parse(file_path, "fasta"):
            sequences.append(str(record.seq).upper())
            ids.append(record.id)
    else:  # text file with one sequence per line
        with open(file_path, 'r') as f:
            for i, line in enumerate(f):
                seq = line.strip().upper()
                if seq:  # Skip empty lines
                    sequences.append(seq)
                    ids.append(f"seq_{i+1}")
    
    return sequences, ids

def setup_dnabert_model(device):
    """Set up DNABERT model."""
    try:
        from transformers import AutoTokenizer, AutoModel
        logger.info("Loading DNABERT model")
        
        # DNABERT-2 is the latest version
        model_name = "zhihan1996/DNABERT-2"
        tokenizer = AutoTokenizer.from_pretrained(model_name)
        model = AutoModel.from_pretrained(model_name)
        
        model = model.to(device)
        model.eval()
        
        return model, tokenizer
    except ImportError:
        logger.error("Transformers module not found. Install with: pip install transformers")
        sys.exit(1)

def setup_nucleotide_transformer_model(device):
    """Set up Nucleotide Transformer model."""
    try:
        from transformers import AutoTokenizer, AutoModelForMaskedLM
        logger.info("Loading Nucleotide Transformer model")
        
        model_name = "InstaDeepAI/nucleotide-transformer-500m-human-ref"
        tokenizer = AutoTokenizer.from_pretrained(model_name)
        model = AutoModelForMaskedLM.from_pretrained(model_name)
        
        model = model.to(device)
        model.eval()
        
        return model, tokenizer
    except ImportError:
        logger.error("Transformers module not found. Install with: pip install transformers")
        sys.exit(1)

def setup_dna_vae_model(input_dim, latent_dim, device):
    """Set up VAE model for DNA embedding generation and anomaly detection."""
    try:
        import torch.nn as nn
        import torch.nn.functional as F
        
        class DNAVAE(nn.Module):
            def __init__(self, input_dim, hidden_dim=512, latent_dim=64):
                super(DNAVAE, self).__init__()
                
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
        
        logger.info(f"Creating DNA VAE model with input dim {input_dim} and latent dim {latent_dim}")
        model = DNAVAE(input_dim, hidden_dim=512, latent_dim=latent_dim)
        model = model.to(device)
        return model
    
    except ImportError:
        logger.error("PyTorch not properly installed. Check your installation.")
        sys.exit(1)

def setup_dna_deep_svdd_model(input_dim, latent_dim, device):
    """Set up Deep SVDD model for DNA anomaly detection."""
    try:
        import torch.nn as nn
        import torch.nn.functional as F
        
        class DNADeepSVDD(nn.Module):
            def __init__(self, input_dim, hidden_dim=512, output_dim=64):
                super(DNADeepSVDD, self).__init__()
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
        
        logger.info(f"Creating DNA DeepSVDD model with input dim {input_dim} and output dim {latent_dim}")
        model = DNADeepSVDD(input_dim, hidden_dim=512, output_dim=latent_dim)
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

def generate_kmers(sequence, k):
    """Generate k-mers from a sequence."""
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

def compute_kmer_frequencies(sequences, k):
    """Compute k-mer frequencies for sequences."""
    # First pass: gather all possible k-mers
    all_kmers = set()
    for seq in sequences:
        kmers = generate_kmers(seq, k)
        all_kmers.update(kmers)
    
    # Create mapping from k-mer to index
    kmer_to_idx = {kmer: i for i, kmer in enumerate(sorted(all_kmers))}
    
    # Second pass: compute frequency vectors
    embeddings = []
    for seq in tqdm(sequences, desc="Computing k-mer frequencies"):
        kmers = generate_kmers(seq, k)
        
        # Count frequencies
        freq = defaultdict(int)
        for kmer in kmers:
            freq[kmer] += 1
        
        # Normalize by sequence length
        for kmer in freq:
            freq[kmer] /= len(kmers)
        
        # Create embedding vector
        vec = np.zeros(len(kmer_to_idx))
        for kmer, count in freq.items():
            if kmer in kmer_to_idx:  # Handle case of rare k-mers not in reference
                vec[kmer_to_idx[kmer]] = count
        
        embeddings.append(vec)
    
    return np.array(embeddings), kmer_to_idx

def process_sequence_with_sliding_window(seq, model, tokenizer, device, window_size, step_size):
    """Process long sequences using sliding window approach."""
    embeddings = []
    
    for i in range(0, len(seq) - window_size + 1, step_size):
        window = seq[i:i+window_size]
        
        inputs = tokenizer(window, return_tensors="pt", truncation=True)
        inputs = {k: v.to(device) for k, v in inputs.items()}
        
        with torch.no_grad():
            outputs = model(**inputs)
            
            # Use [CLS] token embedding or average of all tokens
            if hasattr(outputs, "last_hidden_state"):
                # Take mean of all token embeddings (excluding special tokens)
                embeddings.append(outputs.last_hidden_state[0, 1:-1].mean(dim=0).cpu().numpy())
            else:
                # Fallback for different model architectures
                embeddings.append(outputs[0][0].mean(dim=0).cpu().numpy())
    
    # Aggregate window embeddings (average)
    if embeddings:
        return np.mean(embeddings, axis=0)
    else:
        logger.warning(f"No embeddings generated for sequence of length {len(seq)}")
        return np.zeros(768)  # Default embedding dimension

def generate_dnabert_embeddings(model, tokenizer, sequences, ids, batch_size, device, max_seq_len, window_size, step_size):
    """Generate embeddings using DNABERT model."""
    embeddings = {}
    
    for i in tqdm(range(0, len(sequences), batch_size)):
        batch_ids = ids[i:i+batch_size]
        batch_seqs = sequences[i:i+batch_size]
        
        for j, (id, seq) in enumerate(zip(batch_ids, batch_seqs)):
            if len(seq) <= max_seq_len:
                # For shorter sequences, process directly
                inputs = tokenizer(seq, return_tensors="pt", truncation=True, max_length=max_seq_len)
                inputs = {k: v.to(device) for k, v in inputs.items()}
                
                with torch.no_grad():
                    outputs = model(**inputs)
                    # Average the token embeddings
                    embedding = outputs.last_hidden_state[0].mean(dim=0).cpu().numpy()
                    embeddings[id] = embedding
            else:
                # For longer sequences, use sliding window
                embeddings[id] = process_sequence_with_sliding_window(
                    seq, model, tokenizer, device, window_size, step_size
                )
    
    return embeddings

def generate_nucleotide_transformer_embeddings(model, tokenizer, sequences, ids, batch_size, device, max_seq_len, window_size, step_size):
    """Generate embeddings using Nucleotide Transformer model."""
    embeddings = {}
    
    for i in tqdm(range(0, len(sequences), batch_size)):
        batch_ids = ids[i:i+batch_size]
        batch_seqs = sequences[i:i+batch_size]
        
        for j, (id, seq) in enumerate(zip(batch_ids, batch_seqs)):
            if len(seq) <= max_seq_len:
                # For shorter sequences, process directly
                inputs = tokenizer(seq, return_tensors="pt", truncation=True, max_length=max_seq_len)
                inputs = {k: v.to(device) for k, v in inputs.items()}
                
                with torch.no_grad():
                    outputs = model(**inputs, output_hidden_states=True)
                    # Use the last hidden layer
                    last_hidden_state = outputs.hidden_states[-1]
                    # Average the token embeddings (exclude special tokens)
                    embedding = last_hidden_state[0, 1:-1].mean(dim=0).cpu().numpy()
                    embeddings[id] = embedding
            else:
                # For longer sequences, use sliding window
                embeddings[id] = process_sequence_with_sliding_window(
                    seq, model, tokenizer, device, window_size, step_size
                )
    
    return embeddings

def generate_dna_vae_embeddings(input_features, ids, latent_dim, epochs, learning_rate, batch_size, device):
    """Generate embeddings using VAE model."""
    import torch
    import torch.utils.data as data_utils
    
    # Convert input features to torch dataset
    input_features_tensor = torch.tensor(input_features, dtype=torch.float32)
    dataset = data_utils.TensorDataset(input_features_tensor)
    dataloader = data_utils.DataLoader(dataset, batch_size=batch_size, shuffle=True)
    
    # Create and train VAE model
    input_dim = input_features.shape[1]
    model = setup_dna_vae_model(input_dim, latent_dim, device)
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

def generate_dna_deep_svdd_embeddings(input_features, ids, latent_dim, epochs, learning_rate, batch_size, device):
    """Generate embeddings using Deep SVDD model."""
    import torch
    import torch.utils.data as data_utils
    
    # Convert input features to torch dataset
    input_features_tensor = torch.tensor(input_features, dtype=torch.float32)
    dataset = data_utils.TensorDataset(input_features_tensor)
    dataloader = data_utils.DataLoader(dataset, batch_size=batch_size, shuffle=True)
    
    # Create and train Deep SVDD model
    input_dim = input_features.shape[1]
    model = setup_dna_deep_svdd_model(input_dim, latent_dim, device)
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
    sequences, ids = load_sequences(args.input, args.input_type)
    logger.info(f"Loaded {len(sequences)} DNA sequences")
    
    # Filter sequences that are too short
    filtered_seqs = []
    filtered_ids = []
    for seq, id in zip(sequences, ids):
        if len(seq) > 0:
            # Replace non-ACGT characters with N
            cleaned_seq = ''.join('N' if base not in 'ACGT' else base for base in seq.upper())
            filtered_seqs.append(cleaned_seq)
            filtered_ids.append(id)
        else:
            logger.warning(f"Skipping empty sequence {id}")
    
    if not filtered_seqs:
        logger.error("No valid sequences found after filtering")
        sys.exit(1)
    
    logger.info(f"After filtering: {len(filtered_seqs)} DNA sequences")
    
    # Check if enough sequences for anomaly detection
    if args.detect_anomalies and len(filtered_seqs) < 10:
        logger.warning("Too few sequences for anomaly detection. Disabling.")
        args.detect_anomalies = False
    
    # Generate embeddings based on user choice
    if args.model == 'dnabert':
        model, tokenizer = setup_dnabert_model(args.device)
        embeddings_dict = generate_dnabert_embeddings(
            model, tokenizer, filtered_seqs, filtered_ids, args.batch_size, 
            args.device, args.max_seq_len, args.sliding_window, args.sliding_step
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
                _, anomaly_scores, _ = generate_dna_deep_svdd_embeddings(
                    embeddings_array, ids_list, args.latent_dim, args.epochs, 
                    args.learning_rate, args.batch_size, args.device
                )
                
                # Set anomaly labels based on scores and contamination
                threshold = np.percentile(list(anomaly_scores.values()), 100 * (1 - args.contamination))
                anomaly_labels = {id: score >= threshold for id, score in anomaly_scores.items()}
    
    elif args.model == 'nucleotide_transformer':
        model, tokenizer = setup_nucleotide_transformer_model(args.device)
        embeddings_dict = generate_nucleotide_transformer_embeddings(
            model, tokenizer, filtered_seqs, filtered_ids, args.batch_size, 
            args.device, args.max_seq_len, args.sliding_window, args.sliding_step
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
                _, anomaly_scores, _ = generate_dna_deep_svdd_embeddings(
                    embeddings_array, ids_list, args.latent_dim, args.epochs, 
                    args.learning_rate, args.batch_size, args.device
                )
                
                # Set anomaly labels based on scores and contamination
                threshold = np.percentile(list(anomaly_scores.values()), 100 * (1 - args.contamination))
                anomaly_labels = {id: score >= threshold for id, score in anomaly_scores.items()}
    
    elif args.model == 'k_mer':
        logger.info(f"Using k-mer based embedding with k={args.kmer_size}")
        # Compute k-mer frequencies
        kmer_features, _ = compute_kmer_frequencies(filtered_seqs, args.kmer_size)
        embeddings_dict = {id: emb for id, emb in zip(filtered_ids, kmer_features)}
        
        # Detect anomalies if requested
        if args.detect_anomalies:
            if args.anomaly_method == 'isolation_forest':
                anomaly_labels, anomaly_scores = detect_anomalies_isolation_forest(
                    embeddings_dict, args.contamination
                )
            elif args.anomaly_method == 'vae':
                # Apply VAE to k-mer features
                embeddings_dict, anomaly_scores, _ = generate_dna_vae_embeddings(
                    kmer_features, filtered_ids, args.latent_dim, args.epochs, 
                    args.learning_rate, args.batch_size, args.device
                )
                
                # Set anomaly labels based on scores and contamination
                threshold = np.percentile(list(anomaly_scores.values()), 100 * (1 - args.contamination))
                anomaly_labels = {id: score >= threshold for id, score in anomaly_scores.items()}
            elif args.anomaly_method == 'deep_svdd':
                # Apply DeepSVDD to k-mer features
                embeddings_dict, anomaly_scores, _ = generate_dna_deep_svdd_embeddings(
                    kmer_features, filtered_ids, args.latent_dim, args.epochs, 
                    args.learning_rate, args.batch_size, args.device
                )
                
                # Set anomaly labels based on scores and contamination
                threshold = np.percentile(list(anomaly_scores.values()), 100 * (1 - args.contamination))
                anomaly_labels = {id: score >= threshold for id, score in anomaly_scores.items()}
    
    elif args.model == 'vae':
        logger.info(f"Using VAE for embedding generation with k={args.kmer_size} features")
        # First compute k-mer features as input to VAE
        kmer_features, _ = compute_kmer_frequencies(filtered_seqs, args.kmer_size)
        
        # Generate embeddings using VAE
        embeddings_dict, anomaly_scores, _ = generate_dna_vae_embeddings(
            kmer_features, filtered_ids, args.latent_dim, args.epochs, 
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
            'num_anomalies': sum(1 for val in anomaly_labels.values() if val),
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
    
    logger.info(f"Saving sequence IDs to {args.ids}")
    with open(args.ids, 'w') as f:
        for id in ordered_ids:
            f.write(f"{id}\n")
    
    logger.info(f"Successfully generated embeddings with shape {ordered_embeddings.shape}")

if __name__ == "__main__":
    main()#!/usr/bin/env python
"""
Generate DNA embeddings from FASTA sequences using various embedding models.
Models supported: DNABERT, Nucleotide Transformer, K-mer based
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
from collections import defaultdict
warnings.filterwarnings("ignore")

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('dna_embedding')

def parse_args():
    parser = argparse.ArgumentParser(description='Generate DNA embeddings')
    parser.add_argument('--input', required=True, help='Input DNA sequences file (FASTA or text)')
    parser.add_argument('--input_type', default='text', choices=['fasta', 'text'],
                        help='Input file type (FASTA or line-separated sequences)')
    parser.add_argument('--output', required=True, help='Output NPY file for embeddings')
    parser.add_argument('--ids', required=True, help='Output text file for sequence IDs')
    parser.add_argument('--model', default='dnabert', choices=['dnabert', 'nucleotide_transformer', 'k_mer'],
                        help='Embedding model to use')
    parser.add_argument('--batch_size', type=int, default=8, 
                        help='Batch size for embedding generation')
    parser.add_argument('--max_seq_len', type=int, default=2048, 
                        help='Maximum sequence length')
    parser.add_argument('--kmer_size', type=int, default=6,
                        help='K-mer size for k-mer based embeddings')
    parser.add_argument('--threads', type=int, default=8, 
                        help='Number of CPU threads')
    parser.add_argument('--device', default='cpu', choices=['cpu', 'cuda'],
                        help='Device to use for computation')
    parser.add_argument('--sliding_window', type=int, default=500,
                        help='Sliding window size for long sequences')
    parser.add_argument('--sliding_step', type=int, default=250,
                        help='Step size for sliding window')
    return parser.parse_args()

def load_sequences(file_path, input_type):
    """Load DNA sequences from file."""
    logger.info(f"Loading sequences from {file_path}")
    sequences = []
    ids = []
    
    if input_type == 'fasta':
        for record in SeqIO.parse(file_path, "fasta"):
            sequences.append(str(record.seq).upper())
            ids.append(record.id)
    else:  # text file with one sequence per line
        with open(file_path, 'r') as f:
            for i, line in enumerate(f):
                seq = line.strip().upper()
                if seq:  # Skip empty lines
                    sequences.append(seq)
                    ids.append(f"seq_{i+1}")
    
    return sequences, ids

def setup_dnabert_model(device):
    """Set up DNABERT model."""
    try:
        from transformers import AutoTokenizer, AutoModel
        logger.info("Loading DNABERT model")
        
        # DNABERT-2 is the latest version
        model_name = "zhihan1996/DNABERT-2"
        tokenizer = AutoTokenizer.from_pretrained(model_name)
        model = AutoModel.from_pretrained(model_name)
        
        model = model.to(device)
        model.eval()
        
        return model, tokenizer
    except ImportError:
        logger.error("Transformers module not found. Install with: pip install transformers")
        sys.exit(1)

def setup_nucleotide_transformer_model(device):
    """Set up Nucleotide Transformer model."""
    try:
        from transformers import AutoTokenizer, AutoModelForMaskedLM
        logger.info("Loading Nucleotide Transformer model")
        
        model_name = "InstaDeepAI/nucleotide-transformer-500m-human-ref"
        tokenizer = AutoTokenizer.from_pretrained(model_name)
        model = AutoModelForMaskedLM.from_pretrained(model_name)
        
        model = model.to(device)
        model.eval()
        
        return model, tokenizer
    except ImportError:
        logger.error("Transformers module not found. Install with: pip install transformers")
        sys.exit(1)

def generate_kmers(sequence, k):
    """Generate k-mers from a sequence."""
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

def compute_kmer_frequencies(sequences, k):
    """Compute k-mer frequencies for sequences."""
    # First pass: gather all possible k-mers
    all_kmers = set()
    for seq in sequences:
        kmers = generate_kmers(seq, k)
        all_kmers.update(kmers)
    
    # Create mapping from k-mer to index
    kmer_to_idx = {kmer: i for i, kmer in enumerate(sorted(all_kmers))}
    
    # Second pass: compute frequency vectors
    embeddings = []
    for seq in tqdm(sequences, desc="Computing k-mer frequencies"):
        kmers = generate_kmers(seq, k)
        
        # Count frequencies
        freq = defaultdict(int)
        for kmer in kmers:
            freq[kmer] += 1
        
        # Normalize by sequence length
        for kmer in freq:
            freq[kmer] /= len(kmers)
        
        # Create embedding vector
        vec = np.zeros(len(kmer_to_idx))
        for kmer, count in freq.items():
            if kmer in kmer_to_idx:  # Handle case of rare k-mers not in reference
                vec[kmer_to_idx[kmer]] = count
        
        embeddings.append(vec)
    
    return np.array(embeddings)

def process_sequence_with_sliding_window(seq, model, tokenizer, device, window_size, step_size):
    """Process long sequences using sliding window approach."""
    embeddings = []
    
    for i in range(0, len(seq) - window_size + 1, step_size):
        window = seq[i:i+window_size]
        
        inputs = tokenizer(window, return_tensors="pt", truncation=True)
        inputs = {k: v.to(device) for k, v in inputs.items()}
        
        with torch.no_grad():
            outputs = model(**inputs)
            
            # Use [CLS] token embedding or average of all tokens
            if hasattr(outputs, "last_hidden_state"):
                # Take mean of all token embeddings (excluding special tokens)
                embeddings.append(outputs.last_hidden_state[0, 1:-1].mean(dim=0).cpu().numpy())
            else:
                # Fallback for different model architectures
                embeddings.append(outputs[0][0].mean(dim=0).cpu().numpy())
    
    # Aggregate window embeddings (average)
    if embeddings:
        return np.mean(embeddings, axis=0)
    else:
        logger.warning(f"No embeddings generated for sequence of length {len(seq)}")
        return np.zeros(768)  # Default embedding dimension

def generate_dnabert_embeddings(model, tokenizer, sequences, ids, batch_size, device, max_seq_len, window_size, step_size):
    """Generate embeddings using DNABERT model."""
    embeddings = {}
    
    for i in tqdm(range(0, len(sequences), batch_size)):
        batch_ids = ids[i:i+batch_size]
        batch_seqs = sequences[i:i+batch_size]
        
        for j, (id, seq) in enumerate(zip(batch_ids, batch_seqs)):
            if len(seq) <= max_seq_len:
                # For shorter sequences, process directly
                inputs = tokenizer(seq, return_tensors="pt", truncation=True, max_length=max_seq_len)
                inputs = {k: v.to(device) for k, v in inputs.items()}
                
                with torch.no_grad():
                    outputs = model(**inputs)
                    # Average the token embeddings
                    embedding = outputs.last_hidden_state[0].mean(dim=0).cpu().numpy()
                    embeddings[id] = embedding
            else:
                # For longer sequences, use sliding window
                embeddings[id] = process_sequence_with_sliding_window(
                    seq, model, tokenizer, device, window_size, step_size
                )
    
    return embeddings

def generate_nucleotide_transformer_embeddings(model, tokenizer, sequences, ids, batch_size, device, max_seq_len, window_size, step_size):
    """Generate embeddings using Nucleotide Transformer model."""
    embeddings = {}
    
    for i in tqdm(range(0, len(sequences), batch_size)):
        batch_ids = ids[i:i+batch_size]
        batch_seqs = sequences[i:i+batch_size]
        
        for j, (id, seq) in enumerate(zip(batch_ids, batch_seqs)):
            if len(seq) <= max_seq_len:
                # For shorter sequences, process directly
                inputs = tokenizer(seq, return_tensors="pt", truncation=True, max_length=max_seq_len)
                inputs = {k: v.to(device) for k, v in inputs.items()}
                
                with torch.no_grad():
                    outputs = model(**inputs, output_hidden_states=True)
                    # Use the last hidden layer
                    last_hidden_state = outputs.hidden_states[-1]
                    # Average the token embeddings (exclude special tokens)
                    embedding = last_hidden_state[0, 1:-1].mean(dim=0).cpu().numpy()
                    embeddings[id] = embedding
            else:
                # For longer sequences, use sliding window
                embeddings[id] = process_sequence_with_sliding_window(
                    seq, model, tokenizer, device, window_size, step_size
                )
    
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
    sequences, ids = load_sequences(args.input, args.input_type)
    logger.info(f"Loaded {len(sequences)} DNA sequences")
    
    # Filter sequences that are too short
    filtered_seqs = []
    filtered_ids = []
    for seq, id in zip(sequences, ids):
        if len(seq) > 0:
            # Replace non-ACGT characters with N
            cleaned_seq = ''.join('N' if base not in 'ACGT' else base for base in seq.upper())
            filtered_seqs.append(cleaned_seq)
            filtered_ids.append(id)
        else:
            logger.warning(f"Skipping empty sequence {id}")
    
    if not filtered_seqs:
        logger.error("No valid sequences found after filtering")
        sys.exit(1)
    
    logger.info(f"After filtering: {len(filtered_seqs)} DNA sequences")
    
    # Generate embeddings based on user choice
    if args.model == 'dnabert':
        model, tokenizer = setup_dnabert_model(args.device)
        embeddings_dict = generate_dnabert_embeddings(
            model, tokenizer, filtered_seqs, filtered_ids, args.batch_size, 
            args.device, args.max_seq_len, args.sliding_window, args.sliding_step
        )
    
    elif args.model == 'nucleotide_transformer':
        model, tokenizer = setup_nucleotide_transformer_model(args.device)
        embeddings_dict = generate_nucleotide_transformer_embeddings(
            model, tokenizer, filtered_seqs, filtered_ids, args.batch_size, 
            args.device, args.max_seq_len, args.sliding_window, args.sliding_step
        )
    
    elif args.model == 'k_mer':
        logger.info(f"Using k-mer based embedding with k={args.kmer_size}")
        # For k-mer, we compute directly without batching
        k_mer_embeddings = compute_kmer_frequencies(filtered_seqs, args.kmer_size)
        embeddings_dict = {id: emb for id, emb in zip(filtered_ids, k_mer_embeddings)}
    
    # Convert dictionary to ordered arrays
    ordered_ids = list(embeddings_dict.keys())
    ordered_embeddings = np.array([embeddings_dict[id] for id in ordered_ids])
    
    # Save embeddings and IDs
    logger.info(f"Saving embeddings to {args.output}")
    np.save(args.output, ordered_embeddings)
    
    logger.info(f"Saving sequence IDs to {args.ids}")
    with open(args.ids, 'w') as f:
        for id in ordered_ids:
            f.write(f"{id}\n")
    
    logger.info(f"Successfully generated embeddings with shape {ordered_embeddings.shape}")

if __name__ == "__main__":
    main()

