#!/usr/bin/env python3

"""
Custom implementation of QIIME2's feature-classifier using scikit-learn for taxonomic assignment.
This script provides a standalone version of the Naive Bayes classifier that can be used
outside of the QIIME2 environment if needed.
"""

import argparse
import pandas as pd
import numpy as np
import pickle
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.naive_bayes import MultinomialNB
from sklearn.pipeline import Pipeline
from Bio import SeqIO
import os
import warnings

warnings.filterwarnings('ignore', category=UserWarning)

def parse_args():
    parser = argparse.ArgumentParser(description='Train and use a Naive Bayes classifier for taxonomic assignment')
    subparsers = parser.add_subparsers(dest='command', help='Command to run')

    # Train classifier
    train_parser = subparsers.add_parser('train', help='Train a new classifier')
    train_parser.add_argument('--reference-seqs', required=True, help='Path to reference sequences (FASTA)')
    train_parser.add_argument('--reference-tax', required=True, help='Path to reference taxonomy (TSV)')
    train_parser.add_argument('--output-classifier', required=True, help='Path to output trained classifier')
    train_parser.add_argument('--kmer-length', type=int, default=8, help='k-mer length (default: 8)')

    # Classify sequences
    classify_parser = subparsers.add_parser('classify', help='Classify sequences using a trained classifier')
    classify_parser.add_argument('--classifier', required=True, help='Path to trained classifier')
    classify_parser.add_argument('--input-seqs', required=True, help='Path to input sequences (FASTA)')
    classify_parser.add_argument('--output-tax', required=True, help='Path to output taxonomy assignments')
    classify_parser.add_argument('--confidence-threshold', type=float, default=0.7, help='Confidence threshold (default: 0.7)')

    return parser.parse_args()

def read_fasta(fasta_path):
    """Read sequences from a FASTA file."""
    sequences = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def read_taxonomy(tax_path):
    """Read taxonomy from a TSV file."""
    if tax_path.endswith('.tsv'):
        taxonomy = pd.read_csv(tax_path, sep='\t', index_col=0, header=0)
        if 'Taxon' in taxonomy.columns:
            return taxonomy['Taxon'].to_dict()
        else:
            return taxonomy.iloc[:, 0].to_dict()
    else:
        # Simple format with sequence ID and taxonomy
        taxonomy = {}
        with open(tax_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    seq_id, tax = line.split('\t')
                    taxonomy[seq_id] = tax
        return taxonomy

def train_classifier(reference_seqs, reference_tax, kmer_length=8):
    """Train a Naive Bayes classifier on reference sequences and taxonomy."""
    # Prepare X and y for training
    X = []
    y = []
    
    for seq_id, sequence in reference_seqs.items():
        if seq_id in reference_tax:
            X.append(sequence)
            y.append(reference_tax[seq_id])
    
    # Create k-mer vectorizer
    kmer_vectorizer = CountVectorizer(
        analyzer='char',
        ngram_range=(kmer_length, kmer_length)
    )
    
    # Create and train the classifier
    classifier = Pipeline([
        ('vectorizer', kmer_vectorizer),
        ('classifier', MultinomialNB())
    ])
    
    classifier.fit(X, y)
    
    return classifier

def classify_sequences(classifier, sequences, confidence_threshold=0.7):
    """Classify sequences using the trained classifier."""
    seq_list = list(sequences.values())
    seq_ids = list(sequences.keys())
    
    # Get predictions and probabilities
    predictions = classifier.predict(seq_list)
    probabilities = classifier.predict_proba(seq_list)
    
    # Get confidence scores
    confidence_scores = np.max(probabilities, axis=1)
    
    # Create results dataframe
    results = pd.DataFrame({
        'Feature ID': seq_ids,
        'Taxon': predictions,
        'Confidence': confidence_scores
    })
    
    # Apply confidence threshold
    results.loc[results['Confidence'] < confidence_threshold, 'Taxon'] = 'Unassigned'
    
    return results

def main():
    args = parse_args()
    
    if args.command == 'train':
        # Read reference data
        print("Reading reference sequences...")
        ref_seqs = read_fasta(args.reference_seqs)
        print(f"Read {len(ref_seqs)} reference sequences")
        
        print("Reading reference taxonomy...")
        ref_tax = read_taxonomy(args.reference_tax)
        print(f"Read {len(ref_tax)} reference taxonomies")
        
        # Train classifier
        print("Training classifier...")
        classifier = train_classifier(ref_seqs, ref_tax, args.kmer_length)
        
        # Save classifier
        print(f"Saving classifier to {args.output_classifier}")
        with open(args.output_classifier, 'wb') as f:
            pickle.dump(classifier, f)
        
        print("Classifier training completed")
        
    elif args.command == 'classify':
        # Load classifier
        print(f"Loading classifier from {args.classifier}")
        with open(args.classifier, 'rb') as f:
            classifier = pickle.load(f)
        
        # Read input sequences
        print("Reading input sequences...")
        input_seqs = read_fasta(args.input_seqs)
        print(f"Read {len(input_seqs)} input sequences")
        
        # Classify sequences
        print("Classifying sequences...")
        results = classify_sequences(classifier, input_seqs, args.confidence_threshold)
        
        # Save results
        print(f"Saving taxonomy assignments to {args.output_tax}")
        results.to_csv(args.output_tax, sep='\t', index=False)
        
        print("Classification completed")
    
    else:
        print("Please specify a command: 'train' or 'classify'")

if __name__ == "__main__":
    main()
