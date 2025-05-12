# Metagenomics Pipeline Optimization with Optuna

This document explains the approach to optimizing a metagenomics pipeline using Optuna, focusing on Phase 1 (parameter sensitivity analysis and initial optimization) for the QC and host removal steps.

## Table of Contents

1. [Introduction](#introduction)
2. [Phase 1 Optimization Approach](#phase-1-optimization-approach)
3. [Parameters Being Optimized](#parameters-being-optimized)
4. [Objective Function Design](#objective-function-design)
5. [Code Structure and Workflow](#code-structure-and-workflow)
6. [Running the Optimizer](#running-the-optimizer)
7. [Interpreting the Results](#interpreting-the-results)
8. [Next Steps](#next-steps)

## Introduction

Optimizing the parameters of a complex bioinformatics pipeline can significantly improve results quality. For metagenomic analysis, the early steps of quality control and host removal are critical as they affect all downstream analyses. Using Optuna, a hyperparameter optimization framework, we can systematically explore the parameter space to find optimal settings.

The optimization process follows a phased approach:
- **Phase 1**: Parameter sensitivity analysis and initial optimization
- **Phase 2**: Module-level optimization (to be implemented later)
- **Phase 3**: Integration optimization (to be implemented later)

This document and code focus on Phase 1 for the QC and host removal modules.

## Phase 1 Optimization Approach

Phase 1 consists of two main steps:

### 1. Parameter Sensitivity Analysis

Before full optimization, we test each parameter individually to understand its impact on the pipeline's performance. This helps us:
- Identify which parameters have the greatest influence on results
- Understand the relationship between parameter values and performance
- Prioritize which parameters to focus on in the full optimization

The sensitivity analysis process:
1. Establish a baseline with default parameters
2. For each parameter, test multiple values while keeping others constant
3. Calculate the performance range for each parameter
4. Rank parameters by their impact on performance

### 2. Optuna-Based Optimization

After sensitivity analysis, we use Optuna to find the optimal parameter combination:
1. Define the parameter search space based on the sensitivity analysis
2. Create an objective function that evaluates pipeline performance
3. Use Optuna's Tree-structured Parzen Estimator (TPE) sampler to efficiently explore the parameter space
4. Prune unpromising trials early to save computation
5. Track the best parameter set and performance score

## Parameters Being Optimized

### Fastp (Quality Control)

| Parameter | Type | Range | Description |
|-----------|------|-------|-------------|
| `fastp_qual` | Integer | 15-30 | Quality score threshold for read filtering |
| `fastp_length` | Integer | 40-100 | Minimum read length |
| `fastp_cut_front` | Boolean | True/False | Whether to trim low-quality bases from read start |
| `fastp_cut_tail` | Boolean | True/False | Whether to trim low-quality bases from read end |

### Bowtie2 (Host Removal)

| Parameter | Type | Range | Description |
|-----------|------|-------|-------------|
| `bowtie_mode` | Categorical | very-fast, fast, sensitive, very-sensitive | Alignment sensitivity preset |
| `bowtie_min_score` | Integer | 20-100 | Minimum alignment score for host mapping |

## Objective Function Design

The objective function evaluates pipeline performance using a combination of metrics:

### Primary Metrics

1. **QC Retention Rate**: Percentage of reads that pass quality filtering
   - Higher is better (we want to keep good reads)

2. **Quality Metrics (Q30 Rate)**: Percentage of bases with quality score ≥ 30
   - Higher is better (we want high-quality reads)

3. **Host Removal Accuracy**: How accurately the pipeline removes the expected fraction of host reads
   - Higher is better (closer to expected host fraction)

### Scoring Formula

The overall objective score is calculated as:

```
Score = 0.7 * (QC Retention + Q30 Quality) + 0.3 * Host Removal Accuracy
```

This balances:
- Retaining as many high-quality reads as possible
- Accurately removing host contamination

A higher score indicates better parameter settings.

## Code Structure and Workflow

The optimization code follows this structure:

1. **Configuration and Preparation**
   - Parse command line arguments
   - Load pipeline configuration
   - Prepare a data subset for faster optimization

2. **Sensitivity Analysis**
   - Run pipeline with baseline parameters
   - Test each parameter separately
   - Calculate and rank parameter sensitivity

3. **Optuna Optimization**
   - Define parameter search space
   - Create objective function
   - Run optimization trials
   - Save best parameters

4. **Result Processing**
   - Extract best parameters
   - Generate updated configuration
   - Create visualization of parameter importance

### Key Functions in the Code

- `prepare_subset_data()`: Creates a smaller dataset for faster optimization
- `run_sensitivity_analysis()`: Tests individual parameter impact
- `evaluate_qc_host_removal()`: Evaluates pipeline performance for a parameter set
- `objective_qc_host()`: Optuna objective function
- `run_phase1_optimization()`: Main function coordinating the optimization process

## Running the Optimizer

To run the Phase 1 optimization:

```bash
python parameter_optimization.py --config config.yaml --output optimization_results --trials 50 --subset 0.1 --workers 4 --phase 1 --module qc_host --reference reference_data.json
```

Parameters:
- `--config`: Path to pipeline configuration file
- `--output`: Directory to store optimization results
- `--trials`: Number of Optuna trials to run
- `--subset`: Fraction of data to use (0.0-1.0)
- `--workers`: Number of parallel workers
- `--phase`: Optimization phase (1=sensitivity analysis)
- `--module`: Pipeline module to optimize (qc_host=QC and host removal)
- `--reference`: Optional reference data for evaluation

##