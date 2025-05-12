#!/bin/bash
# Main execution script for MetaVirome Pipeline

set -e

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colorful headers
print_header() {
    echo -e "\n${BLUE}======================${NC}"
    echo -e "${YELLOW}$1${NC}"
    echo -e "${BLUE}======================${NC}\n"
}

# Function to display help message
show_help() {
    echo "MetaVirome Pipeline - A Nextflow workflow for virus detection in metagenomes"
    echo ""
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  -h, --help              Show this help message"
    echo "  -i, --input PATH        Path to input reads (required)"
    echo "  -g, --host PATH         Path to host genome for filtering (required)"
    echo "  -k, --kraken DB         Path to Kraken2 database (required)"
    echo "  -d, --diamond DB        Path to DIAMOND database (required)"
    echo "  -o, --outdir DIR        Output directory (default: results)"
    echo "  -p, --profile PROFILE   Nextflow profile to use (default: standard)"
    echo "  -c, --cpus INT          Maximum CPUs to use (default: 16)"
    echo "  -m, --memory MEM        Maximum memory to use (default: 128.GB)"
    echo "  -s, --single-end        Input is single-end reads (default: paired-end)"
    echo "  -t, --test              Run with test data"
    echo "  -a, --analyze           Run analysis and report generation after pipeline"
    echo ""
    echo "Example:"
    echo "  $0 --input \"data/samples/*_R{1,2}.fastq.gz\" --host data/reference/host.fasta --kraken data/db/kraken2_viral --diamond data/db/viral_proteins"
}

# Default parameters
INPUT=""
HOST_GENOME=""
KRAKEN_DB=""
DIAMOND_DB=""
OUTPUT="results"
PROFILE="standard"
MAX_CPUS="16"
MAX_MEMORY="128.GB"
SINGLE_END="false"
RUN_TEST="false"
RUN_ANALYSIS="false"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -h|--help)
            show_help
            exit 0
            ;;
        -i|--input)
            INPUT="$2"
            shift 2
            ;;
        -g|--host)
            HOST_GENOME="$2"
            shift 2
            ;;
        -k|--kraken)
            KRAKEN_DB="$2"
            shift 2
            ;;
        -d|--diamond)
            DIAMOND_DB="$2"
            shift 2
            ;;
        -o|--outdir)
            OUTPUT="$2"
            shift 2
            ;;
        -p|--profile)
            PROFILE="$2"
            shift 2
            ;;
        -c|--cpus)
            MAX_CPUS="$2"
            shift 2
            ;;
        -m|--memory)
            MAX_MEMORY="$2"
            shift 2
            ;;
        -s|--single-end)
            SINGLE_END="true"
            shift
            ;;
        -t|--test)
            RUN_TEST="true"
            shift
            ;;
        -a|--analyze)
            RUN_ANALYSIS="true"
            shift
            ;;
        *)
            echo -e "${RED}Error: Unknown option $1${NC}"
            show_help
            exit 1
            ;;
    esac
done

# Check if Nextflow is installed
if ! command -v nextflow &> /dev/null; then
    print_header "Nextflow not found, installing..."
    curl -s https://get.nextflow.io | bash
    chmod +x nextflow
    export PATH=$PATH:$PWD
    echo -e "${GREEN}Nextflow installed successfully${NC}"
fi

# Check if we're running the test data
if [[ "$RUN_TEST" == "true" ]]; then
    print_header "Setting up test environment"
    
    # If setup_test_data.sh exists, run it
    if [[ -f setup_test_data.sh ]]; then
        chmod +x setup_test_data.sh
        ./setup_test_data.sh
    else
        echo -e "${RED}Error: Test script setup_test_data.sh not found${NC}"
        exit 1
    fi
    
    # Run with test profile
    print_header "Running pipeline with test data"
    nextflow run metavirome.nf -profile test -resume
    
    # Optionally run analysis
    if [[ "$RUN_ANALYSIS" == "true" ]]; then
        print_header "Running analysis and report generation"
        
        if [[ -f prepare_report_data.sh && -f generate_report.sh ]]; then
            chmod +x prepare_report_data.sh generate_report.sh
            ./generate_report.sh
        else
            echo -e "${RED}Error: Analysis scripts not found${NC}"
            exit 1
        fi
    fi
    
    print_header "Test run complete"
    echo -e "${GREEN}Check results in the test-results directory${NC}"
    exit 0
fi

# If not running test, check required parameters
if [[ "$RUN_TEST" == "false" ]]; then
    if [[ -z "$INPUT" ]]; then
        echo -e "${RED}Error: Input reads path is required${NC}"
        show_help
        exit 1
    fi
    
    if [[ -z "$HOST_GENOME" ]]; then
        echo -e "${RED}Error: Host genome path is required${NC}"
        show_help
        exit 1
    fi
    
    if [[ -z "$KRAKEN_DB" ]]; then
        echo -e "${RED}Error: Kraken2 database path is required${NC}"
        show_help
        exit 1
    fi
    
    if [[ -z "$DIAMOND_DB" ]]; then
        echo -e "${RED}Error: DIAMOND database path is required${NC}"
        show_help
        exit 1
    fi
fi

# Print run information
print_header "MetaVirome Pipeline Configuration"
echo -e "Input reads:       ${YELLOW}$INPUT${NC}"
echo -e "Host genome:       ${YELLOW}$HOST_GENOME${NC}"
echo -e "Kraken2 database:  ${YELLOW}$KRAKEN_DB${NC}"
echo -e "DIAMOND database:  ${YELLOW}$DIAMOND_DB${NC}"
echo -e "Output directory:  ${YELLOW}$OUTPUT${NC}"
echo -e "Profile:           ${YELLOW}$PROFILE${NC}"
echo -e "Max CPUs:          ${YELLOW}$MAX_CPUS${NC}"
echo -e "Max memory:        ${YELLOW}$MAX_MEMORY${NC}"
echo -e "Single-end:        ${YELLOW}$SINGLE_END${NC}"
echo -e "Run analysis:      ${YELLOW}$RUN_ANALYSIS${NC}"

# Confirm to proceed
echo ""
read -p "Proceed with these settings? [y/N] " -n 1 -r
echo ""
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo -e "${RED}Pipeline execution cancelled${NC}"
    exit 0
fi

# Run the pipeline
print_header "Running MetaVirome Pipeline"
nextflow run metavirome.nf \
    --reads "$INPUT" \
    --host_genome "$HOST_GENOME" \
    --kraken2_db "$KRAKEN_DB" \
    --diamond_db "$DIAMOND_DB" \
    --outdir "$OUTPUT" \
    --max_cpus "$MAX_CPUS" \
    --max_memory "$MAX_MEMORY" \
    --single_end "$SINGLE_END" \
    -profile "$PROFILE" \
    -resume

# Check if pipeline completed successfully
if [ $? -ne 0 ]; then
    print_header "Pipeline execution failed"
    echo -e "${RED}Check error messages above for more information${NC}"
    exit 1
fi

print_header "Pipeline execution completed successfully"
echo -e "${GREEN}Results are available in the $OUTPUT directory${NC}"

# Run analysis if requested
if [[ "$RUN_ANALYSIS" == "true" ]]; then
    print_header "Running analysis and report generation"
    
    if [[ -f prepare_report_data.sh && -f generate_report.sh ]]; then
        chmod +x prepare_report_data.sh generate_report.sh
        
        # Modify report paths to use custom output directory if not default
        if [[ "$OUTPUT" != "results" ]]; then
            sed -i "s|results/|$OUTPUT/|g" prepare_report_data.sh
        fi
        
        # Run the report generation
        ./generate_report.sh
        
        if [ $? -eq 0 ]; then
            print_header "Analysis completed successfully"
            echo -e "${GREEN}Final report is available at: final_report/MetaVirome_Analysis_Report.html${NC}"
        else
            print_header "Analysis failed"
            echo -e "${RED}Check error messages above for more information${NC}"
            exit 1
        fi
    else
        echo -e "${RED}Error: Analysis scripts not found${NC}"
        echo -e "${YELLOW}Skipping analysis...${NC}"
    fi
fi

# Generate execution DAG visualization
if command -v dot &> /dev/null; then
    print_header "Generating pipeline visualization"
    nextflow log -f dag | dot -Tpng > pipeline_execution_dag.png
    echo -e "${GREEN}Pipeline visualization saved to: pipeline_execution_dag.png${NC}"
fi

print_header "MetaVirome Pipeline Workflow Complete"
echo -e "${GREEN}Thank you for using the MetaVirome Pipeline!${NC}"
echo -e "${YELLOW}If you found this pipeline useful, please consider citing:${NC}"
echo "  - Kraken2: Wood, D.E., et al. (2019)"
echo "  - DIAMOND: Buchfink, B., et al. (2021)"
echo "  - metaSPAdes: Nurk, S., et al. (2017)"
echo "  - VirSorter2: Guo, J., et al. (2021)"
echo "  - CheckV: Nayfach, S., et al. (2021)"
echo "  - Pharokka: Michniewski, S., et al. (2022)"
echo "  - vConTACT2: Jang, H.B., et al. (2019)"


