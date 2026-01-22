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
    echo "  -P, --parameters_list <list>        List of fixed population mutation rates (comma-separated)."
    echo
    echo "Additional options:"
    echo "  -o, --oriented <1|0>           Indicates if the SFS is oriented (default: 0)."
    echo "  -b, --blocks <num_blocks>      Number of blocks (default: 1)."
    echo "  -u, --upper_bound <value>      Upper bound for the time grid in Ne(0) generations (default: 1)."
    echo "  -l, --lower_bound <value>      Lower bound for the time grid in Ne(0) generations (default: 1e-4)."
    echo "  -e, --epochs <value>          Maximum number of population size changes (default: 5)."
    echo "  -n, --grid_size <value>        Number of time points between lower and upper bound (default: 35)."
    echo "  -t, --troncation <value>       Troncation value (default: 0)."
    echo "  -r, --repeats <value>          Number of repeats for curbe smoothing (default: 0)"
    echo "  -S, --singleton <1|0>          Enable or disable singleton mode (default: 1)."
    echo "      --help                     Display this help message and exit."
    echo
    echo "Example usage:"
    echo "  $0 --sfs my_sfs_file.txt --prefixe_directory ./output --theta_list \"1000,2000,5000\""
    exit 1
}

# Default values
SFS_FILE=""
OUTPUT_DIR=""
ORIENTED=0
NUM_BLOCKS=1
TRONC=0
GENOME_LENGTH=-1
MUTATION_RATE=-1
GENERATION_TIME=-1
UPPER_BOUND=2.
LOWER_BOUND=1e-4
GRID_SIZE=35
EPOCHS=5
REP=0
SING=1
DELTA_FLAG=False;
THETA_LIST=""   # <- LISTE DE THETAS FIXÉS

# Parse command-line arguments with getopt
ARGS=$(getopt -o "s:p:o:b:L:m:g:u:l:e:r:n:S:t:P:d" \
              -l "sfs:,prefixe_directory:,oriented:,blocks:,genome_length:,mutation_rate:,generation_time:,upper_bound:,lower_bound:,epochs:,grid_size:,singleton:,troncation:,parameters_list:,delta_time,help" \
              -- "$@")

if [ $? -ne 0 ]; then
    echo "Error in parsing arguments." >&2
    usage
    exit 1
fi

eval set -- "$ARGS"

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
        -e|--epochs) EPOCHS="$2"; shift 2;;
        -t|--troncation) TRONC="$2"; shift 2;;
        -r|--repeats) REP="$2"; shift 2;;
        -d|--delta_time) DELTA_FLAG=true; shift 1;;   # ✅ flag booléen sans argument
        -n|--grid_size) GRID_SIZE="$2"; shift 2;;
        -S|--singleton) SING="$2"; shift 2;;
        -P|--parameters_list) THETA_LIST="$2"; shift 2;;
        --help) usage; exit 0;;
        --) shift; break;;
        *) echo "Unknown option: $1" >&2; usage; exit 1;;
    esac
done

# Vérification des arguments obligatoires
if [ -z "$SFS_FILE" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: The -s (SFS file) and -p (output_directory) arguments are required."
    usage
    exit 1
fi

# Affichage des options choisies
echo ""
echo "----option list----"
echo ""
echo "SFS file: $SFS_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Oriented: $ORIENTED"
echo "Number of blocks: $NUM_BLOCKS"
echo "Genome length: $GENOME_LENGTH"
echo "Mutation rate: $MUTATION_RATE"
echo "Generation time: $GENERATION_TIME"
echo "Upper bound: $UPPER_BOUND"
echo "Lower bound: $LOWER_BOUND"
echo "Number of epochs: $EPOCHS"
echo "Number of time points: $GRID_SIZE"
echo "Singleton mode: $SING"
echo "Theta list: $THETA_LIST"
echo "Repeats: $REP"
echo

# Exécution du programme C
echo ">>>> Running C program: inference"
START_TIME_C=$(date +%s)
echo 
# Construction de la commande C
C_CMD="./bin/blockbuster_main --sfs $SFS_FILE -p $OUTPUT_DIR -o $ORIENTED -b $NUM_BLOCKS  -u $UPPER_BOUND -l $LOWER_BOUND -e $EPOCHS -n $GRID_SIZE -S $SING -t $TRONC -m $MUTATION_RATE -L $GENOME_LENGTH -g $GENERATION_TIME -r $REP"

# Ajouter -d si le flag est présent
if [ "$DELTA_FLAG" = true ]; then
    C_CMD="$C_CMD -d"
fi

echo
echo "C command line:"
echo "$C_CMD"
echo

# Exécution directe
eval $C_CMD

END_TIME_C=$(date +%s)
EXEC_TIME_C=$((END_TIME_C - START_TIME_C))
echo 
echo "C program execution time: $EXEC_TIME_C seconds"
echo 

# Exécution du programme Python
echo ">>>> Running Python program: plotting"
echo ""
START_TIME_PYTHON=$(date +%s)

python3 parseandplot.py -i "$OUTPUT_DIR/scenarios.txt" -o "$OUTPUT_DIR" \
    -g "$GENERATION_TIME" -or "$ORIENTED"

if [ $? -ne 0 ]; then
    echo "Error: Python program failed to execute."
    exit 1
fi

END_TIME_PYTHON=$(date +%s)
EXEC_TIME_PYTHON=$((END_TIME_PYTHON - START_TIME_PYTHON))
echo "Python program execution time: $EXEC_TIME_PYTHON seconds"
