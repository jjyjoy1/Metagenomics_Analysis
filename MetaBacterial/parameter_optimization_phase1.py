#!/usr/bin/env python3
"""
Metagenomics Pipeline Parameter Optimization - Phase 1
This script performs parameter sensitivity analysis and optimization for the QC and 
host removal steps in the metagenomics pipeline using Optuna.
"""

import os
import sys
import argparse
import subprocess
import pandas as pd
import numpy as np
import optuna
import yaml
import json
import logging
import tempfile
import shutil
from datetime import datetime
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("optimization.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("metagenomic_optimizer")

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Optimize parameters for metagenomic pipeline')
    parser.add_argument('--config', '-c', required=True, 
                        help='Pipeline configuration file (config.yaml)')
    parser.add_argument('--output', '-o', default='optimized_params',
                        help='Directory to store optimization results')
    parser.add_argument('--trials', '-t', type=int, default=50,
                        help='Number of Optuna trials to run')
    parser.add_argument('--subset', '-s', type=float, default=0.1,
                        help='Fraction of data to use for optimization (0.0-1.0)')
    parser.add_argument('--workers', '-w', type=int, default=4,
                        help='Number of parallel workers')
    parser.add_argument('--phase', '-p', type=int, default=1, choices=[1, 2, 3],
                        help='Optimization phase (1=sensitivity, 2=module, 3=integration)')
    parser.add_argument('--module', '-m', default='qc_host',
                        choices=['qc_host', 'taxonomy', 'assembly', 'binning', 'annotation'],
                        help='Pipeline module to optimize')
    parser.add_argument('--reference', '-r', 
                        help='Reference data for evaluation (ground truth)')
    return parser.parse_args()

def load_config(config_file):
    """Load the pipeline configuration file"""
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    return config

def prepare_subset_data(config, subset_fraction, output_dir):
    """Prepare a subset of data for optimization"""
    logger.info(f"Preparing subset of data ({subset_fraction*100:.1f}%) for optimization")
    
    os.makedirs(output_dir, exist_ok=True)
    subset_config = config.copy()
    
    # Create subset for each sample
    subset_samples = {}
    for sample_name, sample_data in config['samples'].items():
        # Create output directories
        sample_subset_dir = os.path.join(output_dir, f"{sample_name}_subset")
        os.makedirs(sample_subset_dir, exist_ok=True)
        
        # Extract subset of reads using seqtk
        r1_subset = os.path.join(sample_subset_dir, f"{sample_name}_R1_subset.fastq.gz")
        r2_subset = os.path.join(sample_subset_dir, f"{sample_name}_R2_subset.fastq.gz")
        
        # Use seqtk to create subset
        seed = 42  # Fixed seed for reproducibility
        subset_size = int(count_reads(sample_data['R1']) * subset_fraction)
        
        if not os.path.exists(r1_subset):
            cmd1 = f"seqtk sample -s{seed} {sample_data['R1']} {subset_size} | gzip > {r1_subset}"
            subprocess.run(cmd1, shell=True, check=True)
            
        if not os.path.exists(r2_subset):
            cmd2 = f"seqtk sample -s{seed} {sample_data['R2']} {subset_size} | gzip > {r2_subset}"
            subprocess.run(cmd2, shell=True, check=True)
        
        # Update config with subset data
        subset_samples[sample_name] = {
            'R1': r1_subset,
            'R2': r2_subset
        }
    
    subset_config['samples'] = subset_samples
    subset_config['outdir'] = os.path.join(output_dir, "subset_results")
    
    # Save subset config
    subset_config_file = os.path.join(output_dir, "subset_config.yaml")
    with open(subset_config_file, 'w') as f:
        yaml.dump(subset_config, f)
    
    return subset_config, subset_config_file

def count_reads(fastq_file):
    """Count number of reads in a FASTQ file"""
    # Assuming 4 lines per read in FASTQ format
    if fastq_file.endswith('.gz'):
        cmd = f"zcat {fastq_file} | wc -l"
    else:
        cmd = f"wc -l {fastq_file}"
    
    result = subprocess.run(cmd, shell=True, check=True, 
                           capture_output=True, text=True)
    line_count = int(result.stdout.strip().split()[0])
    return line_count // 4

def run_snakemake_rule(rule_name, config_file, params=None):
    """Run a specific Snakemake rule with given parameters"""
    cmd = ["snakemake", "--snakefile", "Snakefile", 
           "--configfile", config_file, 
           "--cores", "1", "--use-conda",
           "--nolock", "--rerun-incomplete",
           rule_name]
    
    # Add custom parameters if provided
    if params:
        param_str = ' '.join([f'{k}="{v}"' for k, v in params.items()])
        cmd.extend(["--config", param_str])
    
    logger.debug(f"Running command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        logger.error(f"Snakemake rule {rule_name} failed: {result.stderr}")
        return False
    
    return True

def evaluate_qc_host_removal(config_file, trial_dir, params, reference=None):
    """
    Evaluate the performance of QC and host removal steps
    
    Returns a score based on:
    1. Percentage of non-host reads retained
    2. Percentage of host reads correctly removed
    3. Quality of retained reads
    """
    # Extract parameters
    fastp_params = {
        "qualified_quality_phred": params["fastp_qual"],
        "cut_front": "true" if params["fastp_cut_front"] else "false",
        "cut_tail": "true" if params["fastp_cut_tail"] else "false",
        "length_required": params["fastp_length"]
    }
    
    bowtie_params = {
        "mode": params["bowtie_mode"],
        "min_score": params["bowtie_min_score"]
    }
    
    # Create a temporary config with these parameters
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    # Update config with trial-specific output directory
    config['outdir'] = trial_dir
    os.makedirs(trial_dir, exist_ok=True)
    
    # Write updated config file
    trial_config = os.path.join(trial_dir, "trial_config.yaml")
    with open(trial_config, 'w') as f:
        yaml.dump(config, f)
    
    # Run fastp (QC) step
    fastp_success = run_snakemake_rule("fastp", trial_config, fastp_params)
    if not fastp_success:
        logger.warning("Fastp step failed")
        return -1  # Return a penalty score
    
    # Run bowtie2 (host removal) step
    bowtie_success = run_snakemake_rule("bowtie2_host_removal", trial_config, bowtie_params)
    if not bowtie_success:
        logger.warning("Bowtie2 host removal step failed")
        return -1  # Return a penalty score
    
    # --- Evaluate results ---
    # 1. Calculate read retention rates
    scores = []
    total_host_removal_accuracy = 0
    total_samples = 0
    
    for sample in config['samples']:
        # Check if output files exist
        fastp_dir = os.path.join(trial_dir, "01_fastp")
        hostfree_dir = os.path.join(trial_dir, "02_hostfree")
        
        raw_r1 = config['samples'][sample]['R1']
        trimmed_r1 = os.path.join(fastp_dir, f"{sample}.trimmed.R1.fastq.gz")
        hostfree_r1 = os.path.join(hostfree_dir, f"{sample}.hostfree.R1.fastq.gz")
        
        if not (os.path.exists(trimmed_r1) and os.path.exists(hostfree_r1)):
            logger.warning(f"Missing output files for sample {sample}")
            continue
        
        # Count reads at each stage
        raw_count = count_reads(raw_r1)
        trimmed_count = count_reads(trimmed_r1)
        hostfree_count = count_reads(hostfree_r1)
        
        # Calculate metrics
        qc_retention_rate = trimmed_count / raw_count if raw_count > 0 else 0
        host_removal_rate = (trimmed_count - hostfree_count) / trimmed_count if trimmed_count > 0 else 0
        
        # If we have reference data on expected host content, use it
        if reference and sample in reference:
            expected_host_rate = reference[sample].get('host_fraction', 0.5)  # Default to 50% if not specified
            host_removal_accuracy = 1 - abs(host_removal_rate - expected_host_rate)
            total_host_removal_accuracy += host_removal_accuracy
            total_samples += 1
        
        # Check quality of retained reads from fastp JSON report
        fastp_json = os.path.join(fastp_dir, f"{sample}.fastp.json")
        if os.path.exists(fastp_json):
            with open(fastp_json, 'r') as f:
                fastp_data = json.load(f)
                q30_rate = fastp_data['summary']['after_filtering'].get('q30_rate', 0)
        else:
            q30_rate = 0
        
        # Calculate sample score - balancing retention and quality
        # We want high retention of non-host reads and high removal of host reads
        retention_quality_score = 0.7 * qc_retention_rate + 0.3 * q30_rate
        
        # For now we'll use a simple score
        sample_score = retention_quality_score
        scores.append(sample_score)
    
    # Calculate overall score
    if not scores:
        logger.warning("No valid samples to evaluate")
        return -1
    
    # If we have reference data, factor in host removal accuracy
    if total_samples > 0:
        avg_host_removal_accuracy = total_host_removal_accuracy / total_samples
        overall_score = 0.7 * np.mean(scores) + 0.3 * avg_host_removal_accuracy
    else:
        overall_score = np.mean(scores)
    
    return overall_score

def objective_qc_host(trial, config_file, output_dir, reference=None):
    """Optuna objective function for QC and host removal optimization"""
    # Define parameters to optimize
    params = {
        # Fastp parameters
        "fastp_qual": trial.suggest_int("fastp_qual", 15, 30),
        "fastp_length": trial.suggest_int("fastp_length", 40, 100),
        "fastp_cut_front": trial.suggest_categorical("fastp_cut_front", [True, False]),
        "fastp_cut_tail": trial.suggest_categorical("fastp_cut_tail", [True, False]),
        
        # Bowtie2 parameters
        "bowtie_mode": trial.suggest_categorical("bowtie_mode", 
                                                ["--very-fast", "--fast", "--sensitive", "--very-sensitive"]),
        "bowtie_min_score": trial.suggest_int("bowtie_min_score", 20, 100)
    }
    
    # Create a directory for this trial
    trial_dir = os.path.join(output_dir, f"trial_{trial.number}")
    os.makedirs(trial_dir, exist_ok=True)
    
    # Log the parameters
    logger.info(f"Trial {trial.number} parameters: {params}")
    
    # Run evaluation
    score = evaluate_qc_host_removal(config_file, trial_dir, params, reference)
    
    # Log the result
    logger.info(f"Trial {trial.number} score: {score:.4f}")
    
    return score

def run_sensitivity_analysis(config_file, output_dir, n_trials=10, reference=None):
    """Run parameter sensitivity analysis for Phase 1"""
    logger.info("Running parameter sensitivity analysis (Phase 1)")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    sensitivity_dir = os.path.join(output_dir, "sensitivity_analysis")
    os.makedirs(sensitivity_dir, exist_ok=True)
    
    # Load config
    config = load_config(config_file)
    
    # Define baseline parameters
    baseline_params = {
        "fastp_qual": 20,
        "fastp_length": 50,
        "fastp_cut_front": True,
        "fastp_cut_tail": True,
        "bowtie_mode": "--sensitive",
        "bowtie_min_score": 60
    }
    
    # Run baseline
    logger.info("Running baseline with default parameters")
    baseline_dir = os.path.join(sensitivity_dir, "baseline")
    baseline_score = evaluate_qc_host_removal(config_file, baseline_dir, baseline_params, reference)
    logger.info(f"Baseline score: {baseline_score:.4f}")
    
    # List of parameters to analyze
    params_to_analyze = list(baseline_params.keys())
    sensitivity_results = {}
    
    # Test each parameter independently
    for param in params_to_analyze:
        logger.info(f"Analyzing sensitivity of parameter: {param}")
        param_results = []
        
        # Define range for this parameter
        if param == "fastp_qual":
            values = [15, 20, 25, 30]
        elif param == "fastp_length":
            values = [40, 50, 75, 100]
        elif param == "fastp_cut_front" or param == "fastp_cut_tail":
            values = [True, False]
        elif param == "bowtie_mode":
            values = ["--very-fast", "--fast", "--sensitive", "--very-sensitive"]
        elif param == "bowtie_min_score":
            values = [20, 40, 60, 80, 100]
        else:
            continue
        
        # Test each value
        for value in values:
            # Create a copy of baseline parameters
            test_params = baseline_params.copy()
            test_params[param] = value
            
            # Create directory for this test
            test_dir = os.path.join(sensitivity_dir, f"{param}_{value}")
            
            # Run evaluation
            score = evaluate_qc_host_removal(config_file, test_dir, test_params, reference)
            param_results.append((value, score))
            logger.info(f"  {param} = {value}: score = {score:.4f}")
        
        # Store results
        sensitivity_results[param] = param_results
    
    # Calculate sensitivity scores (how much parameter affects the outcome)
    sensitivity_scores = {}
    for param, results in sensitivity_results.items():
        scores = [score for _, score in results]
        if scores:
            # Calculate range of scores
            score_range = max(scores) - min(scores)
            sensitivity_scores[param] = score_range
    
    # Sort parameters by sensitivity
    sorted_params = sorted(sensitivity_scores.items(), key=lambda x: x[1], reverse=True)
    
    # Save results
    sensitivity_output = {
        "baseline": {
            "parameters": baseline_params,
            "score": baseline_score
        },
        "parameter_sensitivity": {param: score for param, score in sorted_params},
        "detailed_results": {param: results for param, results in sensitivity_results.items()}
    }
    
    # Write results to file
    results_file = os.path.join(output_dir, "sensitivity_results.json")
    with open(results_file, 'w') as f:
        json.dump(sensitivity_output, f, indent=2)
    
    logger.info(f"Sensitivity analysis results saved to {results_file}")
    logger.info("Parameter sensitivity ranking:")
    for param, score in sorted_params:
        logger.info(f"  {param}: {score:.4f}")
    
    return sensitivity_output

def run_phase1_optimization(args):
    """Run Phase 1 optimization (sensitivity analysis and initial optimization)"""
    logger.info(f"Starting Phase 1 optimization for {args.module} module")
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Load config
    config = load_config(args.config)
    
    # Load reference data if provided
    reference = None
    if args.reference:
        with open(args.reference, 'r') as f:
            reference = json.load(f)
    
    # Prepare subset data for optimization
    subset_dir = os.path.join(args.output, "subset_data")
    subset_config, subset_config_file = prepare_subset_data(
        config, args.subset, subset_dir)
    
    # 1. Run sensitivity analysis
    sensitivity_dir = os.path.join(args.output, "sensitivity")
    sensitivity_results = run_sensitivity_analysis(
        subset_config_file, sensitivity_dir, n_trials=10, reference=reference)
    
    # 2. Run Optuna optimization based on sensitive parameters
    optuna_dir = os.path.join(args.output, "optuna_trials")
    os.makedirs(optuna_dir, exist_ok=True)
    
    # Create Optuna study
    logger.info("Starting Optuna optimization")
    study_name = f"phase1_{args.module}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    study_storage = f"sqlite:///{os.path.join(args.output, f'{study_name}.db')}"
    
    if args.module == "qc_host":
        # Create and run study
        study = optuna.create_study(
            study_name=study_name,
            storage=study_storage,
            direction="maximize",
            sampler=optuna.samplers.TPESampler(multivariate=True),
            pruner=optuna.pruners.MedianPruner(n_startup_trials=5, n_warmup_steps=5)
        )
        
        # Run optimization
        study.optimize(
            lambda trial: objective_qc_host(trial, subset_config_file, optuna_dir, reference),
            n_trials=args.trials,
            n_jobs=min(args.workers, args.trials)
        )
        
        # Get best parameters
        best_params = study.best_params
        best_value = study.best_value
        
        logger.info(f"Best parameters: {best_params}")
        logger.info(f"Best score: {best_value:.4f}")
        
        # Save results
        results = {
            "best_parameters": best_params,
            "best_score": best_value,
            "optimization_history": [
                {"trial": t.number, "params": t.params, "score": t.value}
                for t in study.trials
            ],
            "parameter_importance": optuna.importance.get_param_importances(study)
        }
        
        # Save results to file
        results_file = os.path.join(args.output, "optimization_results.json")
        with open(results_file, 'w') as f:
            # Convert values to serializable format
            json_results = {
                "best_parameters": {k: str(v) if isinstance(v, bool) else v 
                                   for k, v in results["best_parameters"].items()},
                "best_score": results["best_score"],
                "optimization_history": [
                    {
                        "trial": t["trial"],
                        "params": {k: str(v) if isinstance(v, bool) else v 
                                  for k, v in t["params"].items()},
                        "score": t["score"]
                    } for t in results["optimization_history"]
                ],
                "parameter_importance": {k: v for k, v in results["parameter_importance"].items()}
            }
            json.dump(json_results, f, indent=2)
        
        logger.info(f"Optimization results saved to {results_file}")
        
        # Generate updated config file with optimized parameters
        optimized_config = config.copy()
        
        # Add optimized parameters to config
        if "fastp" not in optimized_config:
            optimized_config["fastp"] = {}
        
        optimized_config["fastp"]["qualified_quality_phred"] = best_params["fastp_qual"]
        optimized_config["fastp"]["length_required"] = best_params["fastp_length"]
        optimized_config["fastp"]["cut_front"] = best_params["fastp_cut_front"]
        optimized_config["fastp"]["cut_tail"] = best_params["fastp_cut_tail"]
        
        if "bowtie2" not in optimized_config:
            optimized_config["bowtie2"] = {}
        
        optimized_config["bowtie2"]["mode"] = best_params["bowtie_mode"]
        optimized_config["bowtie2"]["min_score"] = best_params["bowtie_min_score"]
        
        # Save optimized config
        optimized_config_file = os.path.join(args.output, "optimized_config.yaml")
        with open(optimized_config_file, 'w') as f:
            yaml.dump(optimized_config, f)
        
        logger.info(f"Optimized configuration saved to {optimized_config_file}")
    
    else:
        logger.warning(f"Optimization for module {args.module} not implemented yet")
    
    return

def main():
    args = parse_args()
    logger.info(f"Starting metagenomics pipeline parameter optimization")
    logger.info(f"Configuration file: {args.config}")
    logger.info(f"Output directory: {args.output}")
    
    # Run the appropriate optimization phase
    if args.phase == 1:
        run_phase1_optimization(args)
    else:
        logger.error(f"Optimization phase {args.phase} not implemented yet")
        sys.exit(1)
    
    logger.info("Optimization completed successfully")

if __name__ == "__main__":
    main()


