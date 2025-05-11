#This is example of cross-sample novelty detection using VAEs model

##Step 1: Prepare your data
First, ensure cross-sample matrices are in the right format:

Tab-separated TSV files
Samples as rows, features as columns
First column is the sample ID

The core matrices we generated earlier should already be in this format:

Pan-genome matrix: {CROSS_SAMPLE_DIR}/pan_genome_matrix.tsv
KO abundance matrix: {CROSS_SAMPLE_DIR}/ko_abundance_matrix.tsv
Taxonomic matrices: {CROSS_SAMPLE_DIR}/genus_abundance_matrix.tsv (or any other level)
Metadata: {CROSS_SAMPLE_DIR}/integrated_metadata.tsv

##Step 2: Run the novelty detection pipeline
Basic usage with all matrices and ensemble model:
'''
bashpython run_novelty_detection.py \
  --pangenome results/cross_sample/pan_genome_matrix.tsv \
  --ko_matrix results/cross_sample/ko_abundance_matrix.tsv \
  --taxonomy results/cross_sample/genus_abundance_matrix.tsv \
  --metadata results/cross_sample/integrated_metadata.tsv \
  --output_dir results/novelty_detection \
  --model ensemble \
  --log_transform \
  --feature_selection \
  --variance_threshold 0.01 \
  --vae_encoder_layers 128,64 \
  --vae_latent_dim 32 \
  --vae_beta 1.0 \
  --vae_kl_anneal \
  --ensemble_weights 0.6,0.4 \
  --top_n 20

To use only VAE model:
bashpython run_novelty_detection.py \
  --ko_matrix results/cross_sample/ko_abundance_matrix.tsv \
  --taxonomy results/cross_sample/species_abundance_matrix.tsv \
  --model vae \
  --vae_encoder_layers 256,128,64 \
  --vae_latent_dim 32 \
  --vae_beta 0.8 \
  --vae_dropout 0.2 \
  --vae_epochs 200 \
  --vae_kl_anneal

To use only Isolation Forest:
bashpython run_novelty_detection.py \
  --pangenome results/cross_sample/pan_genome_matrix.tsv \
  --model isolation_forest \
  --if_n_estimators 200 \
  --if_contamination 0.03
'''

##Step 3: Analyze the results
The script generates several outputs:

anomaly_scores.tsv: Anomaly scores for all samples
top_anomalies.tsv: Top N anomalous samples
annotated_top_anomalies.tsv: Top anomalies with metadata
Various visualizations (score distributions, latent space plots, etc.)

##Next Steps and Customizations

1. Extend the models:
Add more VAE variants (Conditional VAE, β-VAE, etc.)
Implement other scikit-learn models (DBSCAN, One-Class SVM)

2. Add more advanced features:
Implement cross-validation for model evaluation
Add bootstrapping for uncertainty estimation
Implement feature importance analysis

3. Improve visualization:
Add interactive visualization with Plotly
Create a comprehensive dashboard

4. Integration with biological information:
Add taxonomic distance as features
Incorporate functional pathway completeness
Add evolutionary metrics

This implementation provides a solid foundation for novel pathogen detection using both VAEs and an Isolation Forest baseline model. The code is modular, well-documented, and can be extended with additional features as needed.


