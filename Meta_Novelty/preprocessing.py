# preprocessing.py
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.feature_selection import VarianceThreshold

class DataPreprocessor:
    def __init__(self, 
                 scaler_type='standard',
                 feature_selection=True,
                 variance_threshold=0.01,
                 log_transform=True,
                 pca_components=None):
        """
        Initialize data preprocessor for VAE input.
        
        Args:
            scaler_type: 'standard' or 'minmax'
            feature_selection: Whether to perform feature selection
            variance_threshold: Threshold for variance filtering
            log_transform: Whether to apply log transform to count data
            pca_components: Number of PCA components (None to skip PCA)
        """
        self.scaler_type = scaler_type
        self.feature_selection = feature_selection
        self.variance_threshold = variance_threshold
        self.log_transform = log_transform
        self.pca_components = pca_components
        
        # Initialize components
        if scaler_type == 'standard':
            self.scaler = StandardScaler()
        else:
            self.scaler = MinMaxScaler()
            
        self.var_selector = VarianceThreshold(threshold=variance_threshold)
        self.imputer = SimpleImputer(strategy='mean')
        self.pca = None if pca_components is None else PCA(n_components=pca_components)
        
        # Track selected features
        self.selected_features = None
        
    def fit_transform(self, X, matrix_type='count'):
        """
        Fit the preprocessing pipeline and transform the data.
        
        Args:
            X: pandas DataFrame with features
            matrix_type: 'count', 'binary', or 'compositional'
            
        Returns:
            Transformed numpy array
        """
        # Store original feature names
        self.original_features = X.columns.tolist()
        
        # Step 1: Handle different matrix types
        if matrix_type == 'count' and self.log_transform:
            # Log transform for count data (adding pseudo-count to avoid log(0))
            X_processed = np.log1p(X)
        elif matrix_type == 'compositional':
            # Center log-ratio transform for compositional data
            # Add small pseudo-count to avoid log(0)
            X_pseudo = X + 1e-6
            # Compute geometric mean of each row
            g_mean = np.exp(np.mean(np.log(X_pseudo), axis=1))
            # Apply CLR transform
            X_processed = np.log(X_pseudo.div(g_mean, axis=0))
        else:
            X_processed = X.copy()
        
        # Step 2: Impute missing values
        X_imputed = pd.DataFrame(
            self.imputer.fit_transform(X_processed),
            columns=X_processed.columns
        )
        
        # Step 3: Feature selection
        if self.feature_selection:
            self.var_selector.fit(X_imputed)
            support = self.var_selector.get_support()
            self.selected_features = X_imputed.columns[support].tolist()
            X_selected = X_imputed[self.selected_features]
            print(f"Selected {len(self.selected_features)} features out of {len(self.original_features)}")
        else:
            X_selected = X_imputed
            self.selected_features = self.original_features
            
        # Step 4: Scaling
        X_scaled = pd.DataFrame(
            self.scaler.fit_transform(X_selected),
            columns=X_selected.columns
        )
        
        # Step 5: Dimensionality reduction (if requested)
        if self.pca is not None:
            X_final = self.pca.fit_transform(X_scaled)
            print(f"Reduced to {X_final.shape[1]} dimensions with PCA")
            # Save the explained variance
            self.explained_variance_ratio = self.pca.explained_variance_ratio_
            return X_final
        else:
            return X_scaled.values
    
    def transform(self, X):
        """Transform new data using the fitted preprocessor."""
        if self.log_transform:
            X_processed = np.log1p(X)
        else:
            X_processed = X.copy()
            
        X_imputed = pd.DataFrame(
            self.imputer.transform(X_processed),
            columns=X_processed.columns
        )
        
        if self.feature_selection:
            X_selected = X_imputed[self.selected_features]
        else:
            X_selected = X_imputed
            
        X_scaled = pd.DataFrame(
            self.scaler.transform(X_selected),
            columns=X_selected.columns
        )
        
        if self.pca is not None:
            X_final = self.pca.transform(X_scaled)
            return X_final
        else:
            return X_scaled.values
    
    def get_feature_names(self):
        """Return the names of selected features."""
        return self.selected_features


