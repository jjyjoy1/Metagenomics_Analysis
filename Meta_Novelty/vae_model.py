# vae_model.py
import numpy as np
import tensorflow as tf
from tensorflow.keras import layers, Model, optimizers
from tensorflow.keras.callbacks import EarlyStopping, ReduceLROnPlateau
import matplotlib.pyplot as plt
import seaborn as sns
import os

class SamplingLayer(layers.Layer):
    """Sampling layer for VAE to sample from the latent space."""
    
    def call(self, inputs):
        z_mean, z_log_var = inputs
        batch = tf.shape(z_mean)[0]
        dim = tf.shape(z_mean)[1]
        epsilon = tf.keras.backend.random_normal(shape=(batch, dim))
        return z_mean + tf.exp(0.5 * z_log_var) * epsilon

class VAEModel:
    def __init__(self, 
                 input_dim, 
                 encoder_layers=[128, 64], 
                 latent_dim=32, 
                 beta=1.0,
                 dropout_rate=0.1,
                 learning_rate=0.001,
                 activation='relu',
                 use_batch_norm=True):
        """
        Initialize Variational Autoencoder model for novelty detection.
        
        Args:
            input_dim: Dimension of input features
            encoder_layers: List of hidden layer sizes for encoder
            latent_dim: Dimension of latent space
            beta: Weight for KL divergence term
            dropout_rate: Dropout rate
            learning_rate: Learning rate for optimizer
            activation: Activation function
            use_batch_norm: Whether to use batch normalization
        """
        self.input_dim = input_dim
        self.encoder_layers = encoder_layers
        self.latent_dim = latent_dim
        self.beta = beta
        self.dropout_rate = dropout_rate
        self.learning_rate = learning_rate
        self.activation = activation
        self.use_batch_norm = use_batch_norm
        
        # Build the model
        self._build_model()
        
    def _build_model(self):
        """Build the VAE model architecture."""
        # Encoder
        encoder_inputs = layers.Input(shape=(self.input_dim,), name='encoder_input')
        x = encoder_inputs
        
        # Add encoder layers
        for i, units in enumerate(self.encoder_layers):
            x = layers.Dense(units, activation=self.activation, name=f'encoder_dense_{i}')(x)
            
            if self.use_batch_norm:
                x = layers.BatchNormalization(name=f'encoder_bn_{i}')(x)
                
            if self.dropout_rate > 0:
                x = layers.Dropout(self.dropout_rate, name=f'encoder_dropout_{i}')(x)
        
        # Latent space parameters
        z_mean = layers.Dense(self.latent_dim, name='z_mean')(x)
        z_log_var = layers.Dense(self.latent_dim, name='z_log_var')(x)
        
        # Sampling from latent space
        z = SamplingLayer()([z_mean, z_log_var])
        
        # Define encoder model
        self.encoder = Model(encoder_inputs, [z_mean, z_log_var, z], name='encoder')
        
        # Decoder
        decoder_inputs = layers.Input(shape=(self.latent_dim,), name='decoder_input')
        x = decoder_inputs
        
        # Add decoder layers (reverse of encoder)
        for i, units in enumerate(reversed(self.encoder_layers)):
            x = layers.Dense(units, activation=self.activation, name=f'decoder_dense_{i}')(x)
            
            if self.use_batch_norm:
                x = layers.BatchNormalization(name=f'decoder_bn_{i}')(x)
                
            if self.dropout_rate > 0:
                x = layers.Dropout(self.dropout_rate, name=f'decoder_dropout_{i}')(x)
        
        # Output layer
        decoder_outputs = layers.Dense(self.input_dim, name='decoder_output')(x)
        
        # Define decoder model
        self.decoder = Model(decoder_inputs, decoder_outputs, name='decoder')
        
        # Define VAE model
        outputs = self.decoder(self.encoder(encoder_inputs)[2])
        self.vae = Model(encoder_inputs, outputs, name='vae')
        
        # Define custom loss function
        def vae_loss(y_true, y_pred):
            # Reconstruction loss
            reconstruction_loss = tf.reduce_mean(
                tf.keras.losses.mean_squared_error(y_true, y_pred)) * self.input_dim
            
            # KL divergence
            kl_loss = -0.5 * tf.reduce_mean(
                1 + z_log_var - tf.square(z_mean) - tf.exp(z_log_var))
            
            # Total loss
            return reconstruction_loss + self.beta * kl_loss
        
        # Compile the model
        self.vae.compile(
            optimizer=optimizers.Adam(learning_rate=self.learning_rate),
            loss=vae_loss
        )
        
        # Store the losses separately for monitoring
        self.reconstruction_loss_fn = lambda y_true, y_pred: tf.reduce_mean(
            tf.keras.losses.mean_squared_error(y_true, y_pred)) * self.input_dim
        
        self.kl_loss_fn = lambda: -0.5 * tf.reduce_mean(
            1 + z_log_var - tf.square(z_mean) - tf.exp(z_log_var))
    
    def fit(self, X_train, X_val=None, batch_size=32, epochs=100, kl_anneal=True, verbose=1):
        """
        Fit the VAE model.
        
        Args:
            X_train: Training data
            X_val: Validation data (optional)
            batch_size: Batch size
            epochs: Maximum number of epochs
            kl_anneal: Whether to use KL annealing
            verbose: Verbosity level
            
        Returns:
            Training history
        """
        # Set up callbacks
        callbacks = [
            EarlyStopping(monitor='val_loss' if X_val is not None else 'loss', 
                         patience=10, restore_best_weights=True),
            ReduceLROnPlateau(monitor='val_loss' if X_val is not None else 'loss',
                             factor=0.5, patience=5)
        ]
        
        # Implement KL annealing if requested
        if kl_anneal:
            # Custom callback for KL annealing
            class KLAnnealingCallback(tf.keras.callbacks.Callback):
                def __init__(self, beta_start=0, beta_end=1, epochs=epochs):
                    super().__init__()
                    self.beta_start = beta_start
                    self.beta_end = beta_end
                    self.epochs = epochs
                
                def on_epoch_begin(self, epoch, logs=None):
                    # Linear annealing
                    new_beta = self.beta_start + (self.beta_end - self.beta_start) * min(1, epoch / (self.epochs / 2))
                    K = tf.keras.backend
                    K.set_value(self.model.beta, new_beta)
                    print(f"\nEpoch {epoch+1}: beta = {new_beta:.4f}")
            
            # Add the KL annealing callback
            callbacks.append(KLAnnealingCallback())
        
        # Fit the model
        validation_data = (X_val, X_val) if X_val is not None else None
        
        history = self.vae.fit(
            X_train, X_train,
            epochs=epochs,
            batch_size=batch_size,
            shuffle=True,
            validation_data=validation_data,
            callbacks=callbacks,
            verbose=verbose
        )
        
        # Calculate anomaly scores for training data
        self.train_anomaly_scores = self.predict_anomaly_scores(X_train)
        
        # Calculate score statistics
        self.score_mean = np.mean(self.train_anomaly_scores)
        self.score_std = np.std(self.train_anomaly_scores)
        self.score_percentiles = {
            99: np.percentile(self.train_anomaly_scores, 99),
            95: np.percentile(self.train_anomaly_scores, 95),
            90: np.percentile(self.train_anomaly_scores, 90)
        }
        
        return history
    
    def predict_anomaly_scores(self, X, alpha=0.5):
        """
        Predict anomaly scores using the VAE model.
        
        Args:
            X: Input data
            alpha: Weight for reconstruction error vs. KL divergence
                  (alpha=1.0 uses only reconstruction error)
                  
        Returns:
            Anomaly scores (higher = more anomalous)
        """
        # Get latent space representation
        z_mean, z_log_var, _ = self.encoder.predict(X)
        
        # Reconstruct inputs
        X_pred = self.decoder.predict(z_mean)
        
        # Calculate reconstruction error (MSE)
        reconstruction_error = np.mean(np.square(X - X_pred), axis=1)
        
        # Calculate KL divergence
        kl_divergence = -0.5 * np.sum(1 + z_log_var - np.square(z_mean) - np.exp(z_log_var), axis=1)
        
        # Normalize both components
        reconstruction_error_norm = (reconstruction_error - np.mean(reconstruction_error)) / np.std(reconstruction_error)
        kl_divergence_norm = (kl_divergence - np.mean(kl_divergence)) / np.std(kl_divergence)
        
        # Combined anomaly score
        anomaly_scores = alpha * reconstruction_error_norm + (1 - alpha) * kl_divergence_norm
        
        return anomaly_scores
    
    def predict_is_anomaly(self, X, threshold_percentile=95, alpha=0.5):
        """
        Predict if samples are anomalies.
        
        Args:
            X: Input data
            threshold_percentile: Percentile threshold from training scores
            alpha: Weight for reconstruction vs. KL divergence
            
        Returns:
            Boolean array (True = anomaly)
        """
        scores = self.predict_anomaly_scores(X, alpha)
        threshold = self.score_percentiles[threshold_percentile]
        return scores > threshold
    
    def get_latent_representation(self, X):
        """Get the latent space representation for input data."""
        z_mean, _, _ = self.encoder.predict(X)
        return z_mean
    
    def save_model(self, save_dir):
        """Save the VAE model to disk."""
        os.makedirs(save_dir, exist_ok=True)
        
        # Save encoder and decoder models separately
        self.encoder.save(os.path.join(save_dir, 'encoder_model'))
        self.decoder.save(os.path.join(save_dir, 'decoder_model'))
        
        # Save model parameters
        np.savez(os.path.join(save_dir, 'model_params.npz'),
                input_dim=self.input_dim,
                encoder_layers=self.encoder_layers,
                latent_dim=self.latent_dim,
                beta=self.beta,
                dropout_rate=self.dropout_rate,
                learning_rate=self.learning_rate,
                train_score_mean=self.score_mean,
                train_score_std=self.score_std,
                train_score_percentiles=self.score_percentiles)
    
    @classmethod
    def load_model(cls, save_dir):
        """Load the VAE model from disk."""
        # Load model parameters
        params = np.load(os.path.join(save_dir, 'model_params.npz'))
        
        # Create a new instance with the saved parameters
        model = cls(
            input_dim=params['input_dim'].item(),
            encoder_layers=params['encoder_layers'].tolist(),
            latent_dim=params['latent_dim'].item(),
            beta=params['beta'].item(),
            dropout_rate=params['dropout_rate'].item(),
            learning_rate=params['learning_rate'].item()
        )
        
        # Load the trained models
        model.encoder = tf.keras.models.load_model(os.path.join(save_dir, 'encoder_model'))
        model.decoder = tf.keras.models.load_model(os.path.join(save_dir, 'decoder_model'))
        
        # Recreate the full model
        inputs = model.encoder.inputs
        outputs = model.decoder(model.encoder(inputs)[2])
        model.vae = Model(inputs, outputs)
        
        # Load statistics
        model.score_mean = params['train_score_mean'].item()
        model.score_std = params['train_score_std'].item()
        model.score_percentiles = {
            k: params['train_score_percentiles'][i] 
            for i, k in enumerate([99, 95, 90])
        }
        
        return model
    
    def plot_score_distribution(self, scores=None, figsize=(10, 6), bins=50):
        """
        Plot the distribution of anomaly scores.
        
        Args:
            scores: Anomaly scores (if None, uses training scores)
            figsize: Figure size
            bins: Number of histogram bins
        """
        if scores is None:
            scores = self.train_anomaly_scores
            
        plt.figure(figsize=figsize)
        sns.histplot(scores, bins=bins, kde=True)
        
        # Add vertical lines for percentiles
        for percentile, value in self.score_percentiles.items():
            plt.axvline(x=value, color='r', linestyle='--', 
                       label=f'{percentile}th percentile')
            
        plt.title('Distribution of VAE Anomaly Scores')
        plt.xlabel('Anomaly Score (higher = more anomalous)')
        plt.ylabel('Frequency')
        plt.legend()
        plt.tight_layout()
        return plt.gcf()
    
    def plot_latent_space(self, X, labels=None, figsize=(12, 10)):
        """
        Plot samples in 2D latent space.
        
        Args:
            X: Input data
            labels: Labels for coloring points
            figsize: Figure size
        """
        # Get latent space representation
        z_mean = self.get_latent_representation(X)
        
        # If latent dimension is > 2, use PCA to visualize
        if self.latent_dim > 2:
            from sklearn.decomposition import PCA
            pca = PCA(n_components=2)
            z_2d = pca.fit_transform(z_mean)
            title = f'Latent Space (PCA Projection from {self.latent_dim}D)'
        else:
            z_2d = z_mean
            title = '2D Latent Space'
        
        plt.figure(figsize=figsize)
        
        if labels is not None:
            # If labels are continuous, use a colormap
            if np.issubdtype(type(labels[0]), np.number):
                scatter = plt.scatter(z_2d[:, 0], z_2d[:, 1], c=labels, cmap='viridis', alpha=0.7)
                plt.colorbar(scatter, label='Label Value')
            else:
                # For categorical labels
                unique_labels = np.unique(labels)
                for label in unique_labels:
                    indices = np.where(labels == label)[0]
                    plt.scatter(z_2d[indices, 0], z_2d[indices, 1], label=label, alpha=0.7)
                plt.legend()
        else:
            plt.scatter(z_2d[:, 0], z_2d[:, 1], alpha=0.7)
        
        plt.title(title)
        plt.xlabel('Latent Dimension 1')
        plt.ylabel('Latent Dimension 2')
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.tight_layout()
        return plt.gcf()
    
    def plot_reconstruction(self, X, n_samples=5, figsize=(15, 5)):
        """
        Plot original vs reconstructed samples.
        
        Args:
            X: Input data
            n_samples: Number of samples to plot
            figsize: Figure size
        """
        # Only show a subset of features if there are too many
        max_features = 20
        n_features = min(X.shape[1], max_features)
        
        # Choose random samples
        indices = np.random.choice(X.shape[0], n_samples, replace=False)
        X_subset = X[indices, :n_features]
        
        # Reconstruct samples
        X_recon = self.vae.predict(X[indices])[:, :n_features]
        
        # Plot
        plt.figure(figsize=figsize)
        
        for i in range(n_samples):
            # Original
            plt.subplot(2, n_samples, i + 1)
            plt.bar(range(n_features), X_subset[i])
            plt.title(f"Original {i+1}")
            plt.xticks([])
            
            # Reconstruction
            plt.subplot(2, n_samples, n_samples + i + 1)
            plt.bar(range(n_features), X_recon[i])
            plt.title(f"Reconstructed {i+1}")
            plt.xticks([])
        
        plt.tight_layout()
        return plt.gcf()

