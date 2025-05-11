# run_novelty_detection.py
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
import tensorflow as tf
import json
import datetime

# Import our modules
from preprocessing import DataPreprocessor
from isolation_forest import IsolationForestModel
from vae_model import VAEModel
from ensemble_model import EnsembleModel

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Run novelty detection for pathogen discovery')
    
    # Input data arguments
    parser.add_argument('--pangenome', type=str, help='Path to pan-genome matrix TSV')
    parser.add_argument('--ko_matrix', type=str, help='Path to KO abundance matrix TSV')
    parser.add_argument('--taxonomy', type=str, help='Path to taxonomic profile matrix TSV')
    parser.add_argument('--metadata', type=str, help='Path to metadata TSV')
    parser.add_argument('--output_dir', type=str, default='results/novelty_detection',
                       help='Output directory for results')
    
    # Preprocessing arguments
    parser.add_argument('--matrix_type', type=str, default='mixed',
                       choices=['count', 'binary', 'compositional', 'mixed'],
                       help='Type of input matrix')
    parser.add_argument('--log_transform', action='store_true',
                       help='Apply log transform to count data')
    parser.add_argument('--feature_selection', action='store_true',
                       help='Apply feature selection')
    parser.add_argument('--variance_threshold', type=float, default=0.01,
                       help='Variance threshold for feature selection')
    
    # Model arguments
    parser.add_argument('--model', type=str, default='ensemble',
                       choices=['isolation_forest', 'vae', 'ensemble'],
                       help='Which model to use')
    
    # Isolation Forest arguments
    parser.add_argument('--if_n_estimators', type=int, default=100,
                       help='Number of estimators for Isolation Forest')
    parser.add_argument('--if_contamination', type=float, default=0.05,
                       help='Expected contamination for Isolation Forest')
    
    # VAE arguments
    parser.add_argument('--vae_encoder_layers', type=str, default='128,64',
                       help='Comma-separated list of encoder layer sizes for VAE')
    parser.add_argument('--vae_latent_dim', type=int, default=32,
                       help='Latent dimension size for VAE')
    parser.add_argument('--vae_beta', type=float, default=1.0,
                       help='Beta parameter for VAE (KL weight)')
    parser.add_argument('--vae_dropout', type=float, default=0.1,
                       help='Dropout rate for VAE')
    parser.add_argument('--vae_learning_rate', type=float, default=0.001,
                       help='Learning rate for VAE')
    parser.add_argument('--vae_batch_size', type=int, default=32,
                       help='Batch size for VAE training')
    parser.add_argument('--vae_epochs', type=int, default=100,
                       help='Maximum epochs for VAE training')
    parser.add_argument('--vae_kl_anneal', action='store_true',
                       help='Use KL annealing for VAE training')
    
    # Ensemble arguments
    parser.add_argument('--ensemble_weights', type=str, default='0.6,0.4',
                       help='Comma-separated weights for ensemble (VAE,Isolation Forest)')
    
    # Output arguments
    parser.add_argument('--top_n', type=int, default=10,
                       help='Number of top anomalies to report')
    
    return parser.parse_args()

def load_and_preprocess_data(args):
    """Load and preprocess input data."""
    print("Loading input data...")
    
    data_matrices = {}
    sample_ids = None
    
    # Load pan-genome matrix
    if args.pangenome:
        print(f"Loading pan-genome matrix from {args.pangenome}")
        pangenome_df = pd.read_csv(args.pangenome, sep='\t', index_col=0)
        data_matrices['pangenome'] = pangenome_df
        
        # Use these sample IDs as reference
        sample_ids = pangenome_df.index.tolist()
    
    # Load KO matrix
    if args.ko_matrix:
        print(f"Loading KO abundance matrix from {args.ko_matrix}")
        ko_df = pd.read_csv(args.ko_matrix, sep='\t', index_col=0)
        
        # Make sure it has the same samples
        if sample_ids:
            common_samples = set(sample_ids).intersection(ko_df.index)
            if len(common_samples) < len(sample_ids):
                print(f"Warning: KO matrix has {len(common_samples)} out of {len(sample_ids)} samples")
            ko_df = ko_df.loc[sample_ids]
        else:
            sample_ids = ko_df.index.tolist()
            
        data_matrices['ko'] = ko_df
    
    # Load taxonomy matrix
    if args.taxonomy:
        print(f"Loading taxonomy matrix from {args.taxonomy}")
        tax_df = pd.read_csv(args.taxonomy, sep='\t', index_col=0)
        
        # Make sure it has the same samples
        if sample_ids:
            common_samples = set(sample_ids).intersection(tax_df.index)
            if len(common_samples) < len(sample_ids):
                print(f"Warning: Taxonomy matrix has {len(common_samples)} out of {len(sample_ids)} samples")
            tax_df = tax_df.loc[sample_ids]
        else:
            sample_ids = tax_df.index.tolist()
            
        data_matrices['taxonomy'] = tax_df
    
    # Load metadata
    metadata = None
    if args.metadata:
        print(f"Loading metadata from {args.metadata}")
        metadata = pd.read_csv(args.metadata, sep='\t')
        
        # Try to set index to match samples
        id_column = None
        for col in ['sample_id', 'sample', 'id', 'name']:
            if col in metadata.columns:
                id_column = col
                break
                
        if id_column:
            metadata.set_index(id_column, inplace=True)
            
        # Check if metadata covers all samples
        if sample_ids:
            common_samples = set(sample_ids).intersection(metadata.index)
            if len(common_samples) < len(sample_ids):
                print(f"Warning: Metadata covers {len(common_samples)} out of {len(sample_ids)} samples")
                
    # Preprocess each matrix
    preprocessed_data = {}
    preprocessors = {}
    
    print("\nPreprocessing matrices...")
    for name, df in data_matrices.items():
        print(f"Preprocessing {name} matrix...")
        
        # Determine matrix type
        matrix_type = args.matrix_type
        if args.matrix_type == 'mixed':
            if name == 'pangenome':
                matrix_type = 'binary'
            elif name == 'taxonomy':
                matrix_type = 'compositional'
            else:
                matrix_type = 'count'
                
        print(f"  Using {matrix_type} preprocessing for {name} matrix")
        
        # Create and fit preprocessor
        preprocessor = DataPreprocessor(
            scaler_type='standard',
            feature_selection=args.feature_selection,
            variance_threshold=args.variance_threshold,
            log_transform=args.log_transform,
            pca_components=None
        )
        
        # Apply preprocessing
        X_processed = preprocessor.fit_transform(df, matrix_type=matrix_type)
        preprocessed_data[name] = X_processed
        preprocessors[name] = preprocessor
        
        print(f"  {name} matrix processed: {X_processed.shape[0]} samples × {X_processed.shape[1]} features")
    
    # Combine matrices if multiple are provided
    if len(preprocessed_data) > 1:
        print("\nCombining matrices...")
        X_combined = np.concatenate([preprocessed_data[name] for name in preprocessed_data], axis=1)
        print(f"Combined matrix: {X_combined.shape[0]} samples × {X_combined.shape[1]} features")
        
        # Store combined data
        preprocessed_data['combined'] = X_combined
    
    return preprocessed_data, data_matrices, preprocessors, metadata, sample_ids

def train_isolation_forest(X, args):
    """Train Isolation Forest model."""
    print("\nTraining Isolation Forest model...")
    
    model = IsolationForestModel(
        n_estimators=args.if_n_estimators,
        contamination=args.if_contamination,
        random_state=42
    )
    
    model.fit(X)
    print("Isolation Forest model trained")
    
    return model

def train_vae(X, args):
    """Train VAE model."""
    print("\nTraining Variational Autoencoder model...")
    
    # Parse encoder layers
    encoder_layers = [int(x) for x in args.vae_encoder_layers.split(',')]
    
    # Create and train VAE
    vae = VAEModel(
        input_dim=X.shape[1],
        encoder_layers=encoder_layers,
        latent_dim=args.vae_latent_dim,
        beta=args.vae_beta,
        dropout_rate=args.vae_dropout,
        learning_rate=args.vae_learning_rate
    )
    
    # Split data for training/validation
    X_train, X_val = train_test_split(X, test_size=0.2, random_state=42)
    
    # Train the model
    history = vae.fit(
        X_train=X_train,
        X_val=X_val,
        batch_size=args.vae_batch_size,
        epochs=args.vae_epochs,
        kl_anneal=args.vae_kl_anneal
    )
    
    print("VAE model trained")
    
    return vae

def create_ensemble(isolation_forest, vae, args):
    """Create ensemble model."""
    print("\nCreating ensemble model...")
    
    # Parse ensemble weights
    weights_list = [float(x) for x in args.ensemble_weights.split(',')]
    
    # Create ensemble
    models = {
        'vae': vae,
        'isolation_forest': isolation_forest
    }
    
    weights = {
        'vae': weights_list[0],
        'isolation_forest': weights_list[1]
    }
    
    ensemble = EnsembleModel(models, weights)
    print(f"Ensemble created with weights: VAE={weights['vae']:.2f}, IF={weights['isolation_forest']:.2f}")
    
    return ensemble

def evaluate_model(model, X, sample_ids, args, output_dir, model_name):
    """Evaluate model and save results."""
    print(f"\nEvaluating {model_name} model...")
    
    # Create model-specific output directory
    model_dir = os.path.join(output_dir, model_name)
    os.makedirs(model_dir, exist_ok=True)
    
    # Get anomaly scores
    if model_name == 'ensemble':
        scores = model.predict_anomaly_scores(X)
        # We'll use the ensemble scores
        anomaly_scores = scores['ensemble']
    else:
        anomaly_scores = model.predict_anomaly_scores(X)
        
    # Create DataFrame with scores
    results_df = pd.DataFrame({
        'sample_id': sample_ids,
        'anomaly_score': anomaly_scores
    })
    
    # Sort and save results
    results_df = results_df.sort_values('anomaly_score', ascending=False)
    results_path = os.path.join(model_dir, 'anomaly_scores.tsv')
    results_df.to_csv(results_path, sep='\t', index=False)
    print(f"Anomaly scores saved to {results_path}")
    
    # Get top anomalies
    top_anomalies = results_df.head(args.top_n)
    top_path = os.path.join(model_dir, 'top_anomalies.tsv')
    top_anomalies.to_csv(top_path, sep='\t', index=False)
    print(f"Top {args.top_n} anomalies saved to {top_path}")
    
    # Plot score distribution
    if hasattr(model, 'plot_score_distribution'):
        fig = model.plot_score_distribution(anomaly_scores)
        plot_path = os.path.join(model_dir, 'score_distribution.png')
        fig.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Score distribution plot saved to {plot_path}")
    
    # Additional model-specific visualizations
    if model_name == 'vae':
        # Plot latent space
        fig = model.plot_latent_space(X)
        plot_path = os.path.join(model_dir, 'latent_space.png')
        fig.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Latent space plot saved to {plot_path}")
        
        # Plot reconstruction examples
        fig = model.plot_reconstruction(X)
        plot_path = os.path.join(model_dir, 'reconstruction_examples.png')
        fig.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Reconstruction examples saved to {plot_path}")
    
    elif model_name == 'ensemble':
        # Compare models
        fig = model.compare_models(X)
        plot_path = os.path.join(model_dir, 'model_comparison.png')
        fig.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Model comparison plot saved to {plot_path}")
    
    return results_df

def main():
    """Main function to run novelty detection pipeline."""
    # Parse arguments
    args = parse_args()
    
    # Create output directory
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = os.path.join(args.output_dir, f"run_{timestamp}")
    os.makedirs(output_dir, exist_ok=True)
    
    # Save parameters
    with open(os.path.join(output_dir, 'parameters.json'), 'w') as f:
        json.dump(vars(args), f, indent=2)
    
    # Load and preprocess data
    preprocessed_data, raw_data, preprocessors, metadata, sample_ids = load_and_preprocess_data(args)
    
    # Determine which data to use (combined if available, otherwise first available)
    if 'combined' in preprocessed_data:
        X = preprocessed_data['combined']
        data_used = 'combined'
    else:
        data_used = list(preprocessed_data.keys())[0]
        X = preprocessed_data[data_used]
    
    print(f"\nUsing {data_used} data for modeling: {X.shape[0]} samples × {X.shape[1]} features")
    
    # Train models based on user selection
    isolation_forest = None
    vae = None
    ensemble = None
    
    if args.model in ['isolation_forest', 'ensemble']:
        isolation_forest = train_isolation_forest(X, args)
        
        # Save Isolation Forest model
        if isolation_forest is not None:
            os.makedirs(os.path.join(output_dir, 'models'), exist_ok=True)
            isolation_forest.save_model(os.path.join(output_dir, 'models', 'isolation_forest.pkl'))
    
    if args.model in ['vae', 'ensemble']:
        vae = train_vae(X, args)
        
        # Save VAE model
        if vae is not None:
            os.makedirs(os.path.join(output_dir, 'models'), exist_ok=True)
            vae.save_model(os.path.join(output_dir, 'models', 'vae'))
    
    if args.model == 'ensemble':
        ensemble = create_ensemble(isolation_forest, vae, args)
    
    # Evaluate the selected model
    if args.model == 'isolation_forest':
        results = evaluate_model(isolation_forest, X, sample_ids, args, output_dir, 'isolation_forest')
    elif args.model == 'vae':
        results = evaluate_model(vae, X, sample_ids, args, output_dir, 'vae')
    elif args.model == 'ensemble':
        results = evaluate_model(ensemble, X, sample_ids, args, output_dir, 'ensemble')
    
    # If metadata is available, annotate top anomalies with metadata
    if metadata is not None:
        top_anomalies = results.head(args.top_n)
        
        # Try to merge with metadata
        try:
            annotated = pd.merge(
                top_anomalies, 
                metadata, 
                left_on='sample_id', 
                right_index=True,
                how='left'
            )
            
            # Save annotated top anomalies
            annotated_path = os.path.join(output_dir, 'annotated_top_anomalies.tsv')
            annotated.to_csv(annotated_path, sep='\t', index=False)
            print(f"\nAnnotated top anomalies saved to {annotated_path}")
        except Exception as e:
            print(f"\nFailed to annotate top anomalies with metadata: {e}")
    
    print(f"\nNovelty detection completed successfully. Results saved to {output_dir}")

if __name__ == "__main__":
    main()


