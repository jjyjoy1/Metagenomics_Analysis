process VAE_Novelty {
    input:
    path "results/vae_inputs/vae_combined_input.npz"
    path "results/vae_inputs/vae_input_metadata.tsv"

    output:
    path "results/vae_outputs/novelty_scores.tsv"

    script:
    """
    python3 run_vae_novelty.py \
        --input vae_combined_input.npz \
        --metadata vae_input_metadata.tsv \
        --output novelty_scores.tsv
    """
}



