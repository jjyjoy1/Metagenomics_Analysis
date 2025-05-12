# Optuna Parameter Optimization Checklist

This document provides a checklist for implementing Optuna-based parameter optimization with the somatic variant calling pipeline.

## Prerequisites

- [ ] Install Optuna: `pip install optuna`
- [ ] Install visualization dependencies: `pip install plotly matplotlib pandas`
- [ ] Configure a database for storing study results (SQLite, PostgreSQL, etc.)
- [ ] Create a separate configuration file for optimization parameters

## Pipeline Integration Checklist

- [ ] Identify critical parameters to optimize (listed below)
- [ ] Define appropriate parameter ranges for each tool
- [ ] Create objective function that measures pipeline performance
- [ ] Define appropriate metrics (F1 score, precision, recall, etc.)
- [ ] Configure optimization direction (maximize or minimize)
- [ ] Set up appropriate number of trials
- [ ] Implement pruning for inefficient trials
- [ ] Configure parallelization for multiple trials

## Parameters to Optimize

### Read Processing
- [ ] Trimming quality threshold
- [ ] Minimum read length after trimming
- [ ] Adapter stringency

### Alignment (BWA-MEM2)
- [ ] Matching score (`-A`)
- [ ] Mismatch penalty (`-B`) 
- [ ] Gap open penalty (`-O`)
- [ ] Gap extension penalty (`-E`)
- [ ] Seed length (`-k`)
- [ ] Band width (`-w`)

### Duplicate Marking
- [ ] Optical duplicate pixel distance
- [ ] Read name regex for optical duplicates

### BQSR (GATK)
- [ ] Known sites confidence threshold
- [ ] Maximum cycle value

### MuTect2
- [ ] Tumor LOD threshold
- [ ] Normal LOD threshold
- [ ] Downsampling threshold
- [ ] Minimum pruning threshold
- [ ] TLOD (tumor log odds) threshold
- [ ] NLOD (normal log odds) threshold

### VarDict
- [ ] Minimum allele frequency
- [ ] Minimum variant depth
- [ ] Minimum base quality
- [ ] Minimum mapping quality

### FreeBayes
- [ ] Minimum alternate fraction
- [ ] Minimum mapping quality
- [ ] Minimum base quality
- [ ] Minimum alternate count

### Lofreq
- [ ] Minimum coverage
- [ ] Minimum base quality
- [ ] Significance threshold
- [ ] Minimum mapping quality

### Strelka2
- [ ] Maximum indel size
- [ ] Minimum trailing base quality
- [ ] Skip depth filters (on/off)
- [ ] Indel candidate window size

### VarScan2
- [ ] Minimum coverage
- [ ] Minimum variant allele frequency
- [ ] P-value threshold
- [ ] Minimum supporting reads

### DeepVariant
- [ ] Model type selection
- [ ] Number of shards
- [ ] Batch size for inference

### Variant Filtering
- [ ] Min depth threshold
- [ ] Max depth threshold
- [ ] Minimum quality score
- [ ] Strand bias threshold
- [ ] Mapping quality bias threshold
- [ ] Base quality bias threshold
- [ ] Read position bias threshold

### Variant Merging
- [ ] Minimum number of callers in agreement
- [ ] Distance threshold for merging nearby variants

## Performance Evaluation

### Metrics to Consider
- [ ] Precision (minimize false positives)
- [ ] Recall (minimize false negatives)
- [ ] F1 score (balance of precision and recall)
- [ ] Runtime efficiency
- [ ] Memory usage
- [ ] Tumor mutational burden accuracy
- [ ] Concordance with gold standard datasets

### Validation Datasets
- [ ] Synthetic datasets with known variants
- [ ] Gold standard reference samples (e.g., NA12878, SEQC/MAQC)
- [ ] Tumor-normal pairs with orthogonal validation (e.g., PCR-validated variants)
- [ ] Spike-in controls with known variants

## Visualization and Analysis

- [ ] Parameter importance plots
- [ ] Optimization history visualization
- [ ] Parallel coordinate plots
- [ ] Slice plots for key parameters
- [ ] Contour plots for parameter interactions
- [ ] Export of optimal parameter configurations

## Documentation

- [ ] Document baseline performance
- [ ] Document optimization strategy
- [ ] Record performance gains
- [ ] Create configuration template for optimal parameters
- [ ] Provide example optimization script

## Integration with Main Pipeline

- [ ] Create mechanism to apply optimized parameters
- [ ] Add configuration option to use optimized parameters
- [ ] Document how to run pipeline with optimized parameters

## Maintenance

- [ ] Schedule periodic re-optimization
- [ ] Version control for parameter sets
- [ ] Benchmark against new gold standards when available
