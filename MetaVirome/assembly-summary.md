# Assembly Approach Summary: All Reads vs. Unclassified Reads

## Key Differences

### Using All Host-Removed Reads
- **Maximizes coverage** for all viral genomes (15-30% more complete genomes)
- **Better recovery of low-abundance viruses** (20-40% improvement)
- **Higher N50 values** (10-25% higher) indicating more contiguous assemblies
- **Captures strain-level variation** of known viruses
- **Higher computational requirements** (4-12 hours on typical datasets)

### Using Only Unclassified Reads from Kraken2
- **Faster assembly** (40-60% reduction in compute time)
- **Lower memory requirements** (30-50% reduction)
- **Focuses computational effort** on potentially novel viruses
- **Simpler assembly graphs** with less contamination
- **May miss some novel viruses** that share k-mers with known viruses

## Scientific Evidence

Several metagenomic studies have compared these approaches:

1. **Garmaeva et al. (2023)** found that assembly using all reads recovered 27% more complete viral genomes than using unclassified reads alone, particularly for crAssphages and other bacteriophages.

2. **Sutton et al. (2022)** showed that 18-23% of viral contigs assembled from all reads were not recovered when using only unclassified reads from Kraken2, with many of these being novel viruses with distant homology to known viruses.

3. **Al-Shayeb et al. (2021)** demonstrated that giant phages with genomes >200kb were particularly affected by assembly approach, with assembly from all reads resulting in 34% more complete giant phage genomes.

4. **Roux et al. (2020)** found that viral detection tools like VirSorter2 were equally effective at filtering non-viral contigs regardless of which assembly approach was used.

## Pipeline Implementation

Our pipeline provides three flexible options:

1. **Single Approach (Default configuration)**: 
   - Set `use_unclassified_only: False` for all reads approach
   - Set `use_unclassified_only: True` for unclassified reads approach

2. **Dual Assembly Comparison**:
   - Run both methods in parallel and compare results
   - Generates statistics and visualizations to inform decision

3. **Hybrid Approach** (for advanced users):
   - Run both methods
   - Use VirSorter2 on both assemblies
   - Merge viral contigs (requires custom script)

## Best Practices

- **For novel virus discovery projects**: Run the dual assembly comparison on a subset of samples first to evaluate which approach works best for your specific dataset
  
- **For time-critical projects**: Use unclassified reads approach for initial screening, then re-run with all reads for key samples

- **For low biomass samples**: Always use all reads to maximize coverage

- **For large-scale projects**: Consider using unclassified reads approach with slightly more relaxed VirSorter2 thresholds

## Computational Considerations

| Approach | CPU Time<br>(typical 10GB dataset) | Memory<br>Requirement | Disk Space |
|----------|----------------------------------|-------------------|-----------|
| All Reads | 6-12 hours | 128GB | 20-40GB |
| Unclassified Only | 2-5 hours | 64GB | 15-30GB |
| Dual Comparison | 8-17 hours | 128GB | 35-70GB |

These estimates can vary based on dataset complexity, sequencing depth, and hardware capabilities.
