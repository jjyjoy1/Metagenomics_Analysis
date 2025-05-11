# compute_novelty_score.py
# Computes a composite novelty score per contig using multiple features

import pandas as pd
import numpy as np
import argparse


def normalize(series):
    return (series - series.min()) / (series.max() - series.min() + 1e-8)


def compute_score(df, weights):
    df['novelty_score'] = (
        weights['vae_anomaly'] * normalize(df['vae_anomaly']) +
        weights['read_support'] * normalize(df['read_support']) +
        weights['unannotated_fraction'] * normalize(df['unannotated_fraction']) +
        weights['missing_marker_genes'] * normalize(df['missing_marker_genes']) +
        weights['taxonomic_conflict'] * normalize(df['taxonomic_conflict'])
    )
    return df


def main(args):
    df = pd.read_csv(args.input_file, sep='\t')

    weights = {
        'vae_anomaly': args.w1,
        'read_support': args.w2,
        'unannotated_fraction': args.w3,
        'missing_marker_genes': args.w4,
        'taxonomic_conflict': args.w5,
    }

    df = compute_score(df, weights)
    df.sort_values("novelty_score", ascending=False).to_csv(args.output_file, sep='\t', index=False)
    print(f"Saved novelty-ranked contig table to {args.output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', required=True, help='TSV with per-contig metrics')
    parser.add_argument('--output_file', default='ranked_contigs.tsv')
    parser.add_argument('--w1', type=float, default=0.3, help='Weight for VAE anomaly')
    parser.add_argument('--w2', type=float, default=0.2, help='Weight for read support')
    parser.add_argument('--w3', type=float, default=0.2, help='Weight for unannotated fraction')
    parser.add_argument('--w4', type=float, default=0.2, help='Weight for missing marker genes')
    parser.add_argument('--w5', type=float, default=0.1, help='Weight for taxonomic conflict')
    args = parser.parse_args()
    main(args)



