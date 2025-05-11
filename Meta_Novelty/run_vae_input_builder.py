# run_vae_input_builder.py
# Combine k-mer, gene presence, pathway, abundance, and embedding inputs into a single .npz array

import os
import numpy as np
import pandas as pd
import argparse


def load_tsv_or_csv(filepath):
    ext = os.path.splitext(filepath)[1]
    if ext == ".tsv":
        return pd.read_csv(filepath, sep="\t", index_col=0)
    elif ext == ".csv":
        return pd.read_csv(filepath, index_col=0)
    else:
        raise ValueError(f"Unsupported format: {ext}")


def main(args):
    file_map = {}
    for f in os.listdir(args.input_dir):
        for key in ["kmer", "gene_matrix", "pathway", "abundance", "protein_embeddings"]:
            if key in f:
                sid = f.split(f"_{key}")[0]
                file_map.setdefault(sid, {})[key] = os.path.join(args.input_dir, f)

    all_inputs = []
    metadata = []

    for sid, files in file_map.items():
        input_vectors = []
        for key in ["kmer", "gene_matrix", "pathway", "abundance"]:
            if key in files:
                df = load_tsv_or_csv(files[key])
                vec = df.values.flatten()
                input_vectors.append(vec)

        if "protein_embeddings" in files:
            emb = np.load(files["protein_embeddings"])
            input_vectors.append(emb.flatten())

        all_inputs.append(np.concatenate(input_vectors))
        metadata.append({"sample_id": sid, **files})

    all_inputs = np.vstack(all_inputs)
    np.savez(args.output_file, inputs=all_inputs)
    pd.DataFrame(metadata).to_csv(args.metadata_file, sep="\t", index=False)
    print(f"Saved combined input: {args.output_file}\nMetadata: {args.metadata_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", required=True, help="Directory with input .tsv/.csv/.npy files")
    parser.add_argument("--output_file", default="vae_combined_input.npz", help="Output combined .npz file")
    parser.add_argument("--metadata_file", default="vae_input_metadata.tsv", help="Output metadata TSV file")
    args = parser.parse_args()
    main(args)



