#!/bin/bash

# Function to display usage
function usage {
    echo "Usage: $0 -s <sfs> -p <output_directory> [--genome_length <value>] [--mutation_rate <value>] [--generation_time <value>] [options]"
    echo
    echo "Required arguments:"
    echo "  -s, --sfs <sfs>                Path to the Site Frequency Spectrum (SFS) file (required)"
    echo "  -p, --prefixe_directory <dir>  Output directory for the generated results (required)"
    echo
    echo "Optional arguments:"
    echo "  -L, --genome_length <value>    Genome length in number of sites (double, default: -1)"
    echo "  -m, --mutation_rate <value>    Mutation rate per site per generation (double, default: -1)"
    echo "  -g, --generation_time <value>  Generation time in years (double, default: -1)"
    echo
    echo "Additional options:"
    echo "  -o, --oriented <1|0>           Indicates if the SFS is oriented (default: 0). Set to 0 if the SFS is folded."
    echo "  -b, --blocks <num_blocks>      Number of blocks. If b = n, the SFS is divided into training and test sets where the test set contains 1/n sites and the training set contains (n-1)/n sites (default: 1)."
    echo "  -u, --upper_bound <value>      Upper bound for the time grid in Ne(0) generations (default: 1)."
    echo "  -l, --lower_bound <value>      Lower bound for the time grid in Ne(0) generations (default: 1e-4)."
    echo "  -c, --changes <value>          Maximum number of population size changes for which the search is exhaustive. For values between c and 10, a heuristic is applied (default: 5)."
    echo "      --help                     Display this help message and exit."
    echo
    echo "Example usage:"
    echo "  $0 --sfs my_sfs_file.txt --prefixe_directory ./output --genome_length 3000000 --mutation_rate 1e-8 --generation_time 25 --changes 5"
    exit 1
}

# Default values for optional arguments
SFS_FILE=""
OUTPUT_DIR=""
ORIENTED=0
NUM_BLOCKS=1
GENOME_LENGTH=-1
MUTATION_RATE=-1
GENERATION_TIME=-1
UPPER_BOUND=1
LOWER_BOUND=1e-4
CHANGES=5
RECENT="-1"

# Parse command-line arguments with getopt
ARGS=$(getopt -o "s:p:o:b:L:m:g:u:l:c:h:r:" -l "sfs:,prefixe_directory:,oriented:,blocks:,genome_length:,mutation_rate:,generation_time:,upper_bound:,lower_bound:,changes:,help" -- "$@")
if [ $? -ne 0 ]; then
    usage
fi

eval set -- "$ARGS"

# Extract options and their arguments
while true; do
    case "$1" in
        -s|--sfs) SFS_FILE="$2"; shift 2;;
        -p|--prefixe_directory) OUTPUT_DIR="$2"; shift 2;;
        -o|--oriented) ORIENTED="$2"; shift 2;;
        -b|--blocks) NUM_BLOCKS="$2"; shift 2;;
        -L|--genome_length) GENOME_LENGTH="$2"; shift 2;;
        -m|--mutation_rate) MUTATION_RATE="$2"; shift 2;;
        -g|--generation_time) GENERATION_TIME="$2"; shift 2;;
        -u|--upper_bound) UPPER_BOUND="$2"; shift 2;;
        -l|--lower_bound) LOWER_BOUND="$2"; shift 2;;
        -c|--changes) CHANGES="$2"; shift 2;;
        -r|--recent) RECENT="$2"; shift 2;;
        --help) usage; shift;;
        --) shift; break;;
        *) usage; break;;
    esac
done

# Check for required arguments
if [ -z "$SFS_FILE" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: The -s (SFS file) and -p (output_directory) arguments are required."
    usage
fi

# Automatically set ORIENTED to 0 if the SFS file is folded
if [[ "$SFS_FILE" == *"folded"* ]]; then
    ORIENTED=0
    echo "The SFS file is folded. The '-o' option is automatically set to 0."
fi

# Display the chosen options
echo "SFS file: $SFS_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Oriented: $ORIENTED"
echo "Number of blocks: $NUM_BLOCKS"
echo "Genome length: $GENOME_LENGTH"
echo "Mutation rate: $MUTATION_RATE"
echo "Generation time: $GENERATION_TIME"
echo "Upper bound: $UPPER_BOUND"
echo "Lower bound: $LOWER_BOUND"
echo "Number of changes: $CHANGES"

# Measure execution time for C program
echo "Running C program..."
START_TIME_C=$(date +%s)
echo --sfs "$SFS_FILE" -p "$OUTPUT_DIR" -o "$ORIENTED" -b "$NUM_BLOCKS" -u "$UPPER_BOUND" -l "$LOWER_BOUND" -c "$CHANGES" -r "$RECENT"
./bin/blockbuster_main --sfs "$SFS_FILE" -p "$OUTPUT_DIR" -o "$ORIENTED" -b "$NUM_BLOCKS" -u "$UPPER_BOUND" -l "$LOWER_BOUND" -c "$CHANGES" -r "$RECENT"

if [ $? -ne 0 ]; then
    echo "Error: C program failed to execute. Exiting."
    exit 1
fi

END_TIME_C=$(date +%s)
EXEC_TIME_C=$((END_TIME_C - START_TIME_C))
echo "C program execution time: $EXEC_TIME_C seconds"

# Measure execution time for Python program
echo "Running Python program..."
START_TIME_PYTHON=$(date +%s)
python3 parseandplot.py -i "$OUTPUT_DIR/scenarios.txt" -o "$OUTPUT_DIR" -l "$GENOME_LENGTH" -m "$MUTATION_RATE" -g "$GENERATION_TIME"

if [ $? -ne 0 ]; then
    echo "Error: Python program failed to execute."
    exit 1
fi

END_TIME_PYTHON=$(date +%s)
EXEC_TIME_PYTHON=$((END_TIME_PYTHON - START_TIME_PYTHON))
echo "Python program execution time: $EXEC_TIME_PYTHON seconds"
