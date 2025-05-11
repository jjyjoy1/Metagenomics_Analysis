# isolation_forest.py
import numpy as np
import pandas as pd
from sklearn.ensemble import IsolationForest
import joblib
import matplotlib.pyplot as plt
import seaborn as sns

class IsolationForestModel:
    def __init__(self, 
                 n_estimators=100, 
                 contamination=0.05, 
                 random_state=42, 
                 max_features=1.0):
        """
        Initialize Isolation Forest model for novelty detection.
        
        Args:
            n_estimators: Number of trees in the forest
            contamination: Expected proportion of anomalies
            random_state: Random seed
            max_features: Maximum features to consider
        """
        self.n_estimators = n_estimators
        self.contamination = contamination
        self.random_state = random_state
        self.max_features = max_features
        
        self.model = IsolationForest(
            n_estimators=n_estimators,
            contamination=contamination,
            random_state=random_state,
            max_features=max_features,
            n_jobs=-1
        )
    
    def fit(self, X):
        """
        Fit the Isolation Forest model.
        
        Args:
            X: Training data
        """
        self.model.fit(X)
        
        # Calculate anomaly scores for training data
        # Note: Isolation Forest returns negative scores, with more negative = more anomalous
        # We convert to a positive score where higher = more anomalous
        raw_scores = -self.model.decision_function(X)
        self.train_scores = raw_scores
        
        # Calculate statistics for later use
        self.score_mean = np.mean(raw_scores)
        self.score_std = np.std(raw_scores)
        self.score_percentiles = {
            99: np.percentile(raw_scores, 99),
            95: np.percentile(raw_scores, 95),
            90: np.percentile(raw_scores, 90)
        }
        
        return self
    
    def predict_anomaly_scores(self, X):
        """
        Predict anomaly scores for new samples.
        
        Args:
            X: Test data
            
        Returns:
            Anomaly scores (higher = more anomalous)
        """
        return -self.model.decision_function(X)
    
    def predict_is_anomaly(self, X, threshold_percentile=95):
        """
        Predict if samples are anomalies based on percentile threshold.
        
        Args:
            X: Test data
            threshold_percentile: Percentile for threshold (95 = top 5% are anomalies)
            
        Returns:
            Boolean array where True = anomaly
        """
        scores = self.predict_anomaly_scores(X)
        threshold = self.score_percentiles[threshold_percentile]
        return scores > threshold
    
    def get_feature_importances(self, X, feature_names=None):
        """
        Estimate feature importance for anomaly detection (experimental).
        
        This is a heuristic approach since Isolation Forest doesn't provide
        direct feature importances for anomaly detection.
        
        Args:
            X: Input data
            feature_names: Names of features
            
        Returns:
            DataFrame with feature importances
        """
        if feature_names is None:
            feature_names = [f"feature_{i}" for i in range(X.shape[1])]
            
        # We'll use permutation-based approach
        base_scores = self.predict_anomaly_scores(X)
        importances = []
        
        for i in range(X.shape[1]):
            X_permuted = X.copy()
            X_permuted[:, i] = np.random.permutation(X_permuted[:, i])
            permuted_scores = self.predict_anomaly_scores(X_permuted)
            
            # Importance is how much the score changes when we permute a feature
            importance = np.mean(np.abs(permuted_scores - base_scores))
            importances.append(importance)
            
        return pd.DataFrame({
            'feature': feature_names,
            'importance': importances
        }).sort_values('importance', ascending=False)
    
    def save_model(self, filepath):
        """Save the model to disk."""
        joblib.dump(self, filepath)
    
    @classmethod
    def load_model(cls, filepath):
        """Load the model from disk."""
        return joblib.load(filepath)
    
    def plot_score_distribution(self, scores=None, figsize=(10, 6), bins=50):
        """
        Plot the distribution of anomaly scores.
        
        Args:
            scores: Scores to plot (if None, uses training scores)
            figsize: Figure size
            bins: Number of histogram bins
        """
        if scores is None:
            scores = self.train_scores
            
        plt.figure(figsize=figsize)
        sns.histplot(scores, bins=bins, kde=True)
        
        # Add vertical lines for percentiles
        for percentile, value in self.score_percentiles.items():
            plt.axvline(x=value, color='r', linestyle='--', 
                       label=f'{percentile}th percentile')
            
        plt.title('Distribution of Anomaly Scores')
        plt.xlabel('Anomaly Score (higher = more anomalous)')
        plt.ylabel('Frequency')
        plt.legend()
        plt.tight_layout()
        return plt.gcf()


