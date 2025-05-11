# ensemble_model.py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc, precision_recall_curve

class EnsembleModel:
    def __init__(self, models, weights=None):
        """
        Create an ensemble of novelty detection models.
        
        Args:
            models: Dictionary of {'model_name': model_object}
            weights: Dictionary of {'model_name': weight}
                     If None, equal weights are used
        """
        self.models = models
        
        if weights is None:
            # Equal weights
            self.weights = {name: 1.0 / len(models) for name in models}
        else:
            # Normalize weights to sum to 1
            total = sum(weights.values())
            self.weights = {name: weight / total for name, weight in weights.items()}
            
    def predict_anomaly_scores(self, X):
        """
        Get weighted ensemble anomaly scores.
        
        Args:
            X: Input data
            
        Returns:
            Dictionary of {'model_name': scores} and 'ensemble' for combined scores
        """
        scores = {}
        
        # Get scores from each model
        for name, model in self.models.items():
            scores[name] = model.predict_anomaly_scores(X)
            
        # Normalize each model's scores
        normalized_scores = {}
        for name, model_scores in scores.items():
            # Z-score normalization
            normalized_scores[name] = (model_scores - np.mean(model_scores)) / np.std(model_scores)
            
        # Compute weighted ensemble score
        ensemble_scores = np.zeros(X.shape[0])
        for name, model_scores in normalized_scores.items():
            ensemble_scores += self.weights[name] * model_scores
            
        # Add ensemble scores to the dictionary
        scores['ensemble'] = ensemble_scores
        
        return scores
    
    def predict_is_anomaly(self, X, threshold_percentile=95, model_name='ensemble'):
        """
        Predict anomalies using the specified model.
        
        Args:
            X: Input data
            threshold_percentile: Percentile threshold
            model_name: Which model to use ('ensemble' or a specific model name)
            
        Returns:
            Boolean array (True = anomaly)
        """
        scores = self.predict_anomaly_scores(X)
        model_scores = scores[model_name]
        
        # Calculate threshold based on percentile
        threshold = np.percentile(model_scores, threshold_percentile)
        
        return model_scores > threshold
    

    def get_top_anomalies(self, X, sample_ids=None, top_n=10, model_name='ensemble'):
        """
        Get the top anomalies according to the specified model.
        
        Args:
            X: Input data
            sample_ids: Sample identifiers (if None, uses indices)
            top_n: Number of top anomalies to return
            model_name: Which model to use ('ensemble' or a specific model name)
            
        Returns:
            DataFrame with top anomalies and their scores
        """
        scores = self.predict_anomaly_scores(X)
        model_scores = scores[model_name]
        
        # Create sample IDs if not provided
        if sample_ids is None:
            sample_ids = [f"Sample_{i}" for i in range(len(model_scores))]
            
        # Create DataFrame with scores
        results = pd.DataFrame({
            'sample_id': sample_ids,
            'anomaly_score': model_scores
        })
        
        # Sort by score (descending) and return top N
        return results.sort_values('anomaly_score', ascending=False).head(top_n)
    
    def compare_models(self, X, figsize=(12, 8)):
        """
        Compare anomaly scores from different models.
        
        Args:
            X: Input data
            figsize: Figure size
            
        Returns:
            Matplotlib figure
        """
        scores = self.predict_anomaly_scores(X)
        
        # Create a DataFrame with all scores
        df = pd.DataFrame({name: scores[name] for name in scores})
        
        # Compute correlation matrix
        corr_matrix = df.corr()
        
        # Plot correlation heatmap
        plt.figure(figsize=figsize)
        
        # Subplot 1: Correlation heatmap
        plt.subplot(1, 2, 1)
        sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', vmin=-1, vmax=1)
        plt.title('Model Score Correlation')
        
        # Subplot 2: Score distributions
        plt.subplot(1, 2, 2)
        for name in scores:
            sns.kdeplot(scores[name], label=name)
        plt.title('Score Distributions')
        plt.xlabel('Anomaly Score (Normalized)')
        plt.ylabel('Density')
        plt.legend()
        
        plt.tight_layout()
        return plt.gcf()
    
    def plot_roc_curves(self, X, y_true, figsize=(10, 8)):
        """
        Plot ROC curves for each model.
        
        Args:
            X: Input data
            y_true: Ground truth labels (1 for anomaly, 0 for normal)
            figsize: Figure size
            
        Returns:
            Matplotlib figure
        """
        scores = self.predict_anomaly_scores(X)
        
        plt.figure(figsize=figsize)
        
        # Plot ROC curve for each model
        for name, model_scores in scores.items():
            fpr, tpr, _ = roc_curve(y_true, model_scores)
            roc_auc = auc(fpr, tpr)
            plt.plot(fpr, tpr, lw=2, label=f'{name} (AUC = {roc_auc:.3f})')
            
        # Add diagonal line
        plt.plot([0, 1], [0, 1], 'k--', lw=2)
        
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver Operating Characteristic (ROC) Curves')
        plt.legend(loc="lower right")
        plt.grid(True, linestyle='--', alpha=0.7)
        
        return plt.gcf()
    
    def plot_precision_recall_curves(self, X, y_true, figsize=(10, 8)):
        """
        Plot precision-recall curves for each model.
        
        Args:
            X: Input data
            y_true: Ground truth labels (1 for anomaly, 0 for normal)
            figsize: Figure size
            
        Returns:
            Matplotlib figure
        """
        scores = self.predict_anomaly_scores(X)
        
        plt.figure(figsize=figsize)
        
        # Plot precision-recall curve for each model
        for name, model_scores in scores.items():
            precision, recall, _ = precision_recall_curve(y_true, model_scores)
            pr_auc = auc(recall, precision)
            plt.plot(recall, precision, lw=2, label=f'{name} (AUC = {pr_auc:.3f})')
            
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.title('Precision-Recall Curves')
        plt.legend(loc="best")
        plt.grid(True, linestyle='--', alpha=0.7)
        
        return plt.gcf()



