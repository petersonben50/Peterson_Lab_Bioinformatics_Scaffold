#!/bin/bash

#==============================================================================#
#                 PLBS fastp Filtering Executable Wrapper Script               #
#==============================================================================#
# This script serves as a standardized wrapper for running fastp QC/trimming
# within an Apptainer container, adhering to PLBS conventions.
# This read filtering is specifically designed for metagenomic reads that
# will then undergo metagenomic assembly.
# It is designed to be called by Snakemake or directly for testing.

#------------------------------------------------------------------------------#
#                           1. Argument Parsing                                #
#------------------------------------------------------------------------------#
# Parse command-line arguments for customizable parameters.
# This makes the script flexible and reusable.

# Initialize variables with default or empty values
FASTP_CONTAINER=""
METAGENOME_ID=""
NUM_THREADS=""
CONDA_HOME="/opt/conda" # Default Conda path inside container

# Input and Output File Paths
INPUT_READ1=""
INPUT_READ2=""
TRIMMED_R1=""
TRIMMED_R2=""
TRIMMED_UNPAIRED1=""
TRIMMED_UNPAIRED2=""
TRIMMED_MERGED=""
FAILED_READS=""
FASTP_JSON_REPORT=""
FASTP_HTML_REPORT=""

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --fastp_container) FASTP_CONTAINER="$2"; shift ;;
        --metagenome_id) METAGENOME_ID="$2"; shift ;;
        --num_threads) NUM_THREADS="$2"; shift ;;
        --input_read1) INPUT_READ1="$2"; shift ;;
        --input_read2) INPUT_READ2="$2"; shift ;;
        --output_trimmed_r1) TRIMMED_R1="$2"; shift ;;
        --output_trimmed_r2) TRIMMED_R2="$2"; shift ;;
        --output_unpaired1) TRIMMED_UNPAIRED1="$2"; shift ;;
        --output_unpaired2) TRIMMED_UNPAIRED2="$2"; shift ;;
        --output_merged) TRIMMED_MERGED="$2"; shift ;;
        --output_failed) FAILED_READS="$2"; shift ;;
        --output_json) FASTP_JSON_REPORT="$2"; shift ;;
        --output_html) FASTP_HTML_REPORT="$2"; shift ;;
        --conda_home) CONDA_HOME="$2"; shift ;; # Override if different path
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Validate required arguments
if [ -z "$FASTP_CONTAINER" ] || [ -z "$METAGENOME_ID" ] || [ -z "$NUM_THREADS" ] || \
   [ -z "$INPUT_READ1" ] || [ -z "$INPUT_READ2" ] || \
   [ -z "$TRIMMED_R1" ] || [ -z "$TRIMMED_R2" ] || \
   [ -z "$FASTP_JSON_REPORT" ] || [ -z "$FASTP_HTML_REPORT" ]; then
    echo "Error: Missing required arguments."
    echo "Usage: $0 --fastp_container <path> --metagenome_id <id> --num_threads <int> \\"
    echo "          --input_read1 <path> --input_read2 <path> \\"
    echo "          --output_trimmed_r1 <path> --output_trimmed_r2 <path> \\"
    echo "          --output_json <path> --output_html <path> \\"
    echo "          [--output_unpaired1 <path>] [--output_unpaired2 <path>] \\"
    echo "          [--output_merged <path>] [--output_failed <path>] \\"
    echo "          [--conda_home <path>]"
    exit 1
fi

#------------------------------------------------------------------------------#
#                           2. PLBS Path Construction & Validation             #
#------------------------------------------------------------------------------#
# The input directory for input reads.
# Derived from the first input read file to ensure consistency.
INPUT_DIR="$(dirname "${INPUT_READ1}")"

# Validate input files exist
if [ ! -f "${INPUT_READ1}" ]; then
    echo "Error: Input read 1 file does not exist: ${INPUT_READ1}"
    exit 1
fi
if [ ! -f "${INPUT_READ2}" ]; then
    echo "Error: Input read 2 file does not exist: ${INPUT_READ2}"
    exit 1
fi

# The output directory for this sample's fastp results.
# Derived from one of the output files to ensure consistency.
OUTPUT_SAMPLE_DIR="$(dirname "${TRIMMED_R1}")"

# Ensure output directory exists on the host
mkdir -p "${OUTPUT_SAMPLE_DIR}" || { echo "Error: Could not create output directory ${OUTPUT_SAMPLE_DIR}"; exit 1; }


#------------------------------------------------------------------------------#
#                           3. Logging & Pre-execution Checks                  #
#------------------------------------------------------------------------------#
echo "----------------------------------------------------------------"
echo "PLBS fastp QC/Trimming Job Started: $(date)"
echo "Metagenome ID: ${METAGENOME_ID}"
echo "Input R1: ${INPUT_READ1}"
echo "Input R2: ${INPUT_READ2}"
echo "Output Directory: ${OUTPUT_SAMPLE_DIR}"
echo "Fastp Container: ${FASTP_CONTAINER}"
echo "Threads Requested: ${NUM_THREADS}"
echo "----------------------------------------------------------------"

#------------------------------------------------------------------------------#
#                       4. Apptainer Execution of fastp                        #
#------------------------------------------------------------------------------#
# Construct the fastp command string for execution inside the container.
# All paths are relative to the container's bind mounts (/data/input_dir, /data/output_dir).
# Note: Double quotes are escaped with backslashes for the inner bash -c command.
# Also, we now bind the specific input read files, not just their directory.
FASTP_COMMAND_STRING="fastp"
FASTP_COMMAND_STRING+=" --in1 \"/data/input_dir/$(basename "${INPUT_READ1}")\""
FASTP_COMMAND_STRING+=" --in2 \"/data/input_dir/$(basename "${INPUT_READ2}")\""
FASTP_COMMAND_STRING+=" --out1 \"/data/output_dir/$(basename "${TRIMMED_R1}")\""
FASTP_COMMAND_STRING+=" --out2 \"/data/output_dir/$(basename "${TRIMMED_R2}")\""
if [ -n "${TRIMMED_UNPAIRED1}" ]; then
  FASTP_COMMAND_STRING+=" --unpaired1 \"/data/output_dir/$(basename "${TRIMMED_UNPAIRED1}")\""
fi
if [ -n "${TRIMMED_UNPAIRED2}" ]; then
  FASTP_COMMAND_STRING+=" --unpaired2 \"/data/output_dir/$(basename "${TRIMMED_UNPAIRED2}")\""
fi
if [ -n "${TRIMMED_MERGED}" ]; then
  FASTP_COMMAND_STRING+=" --merge --merged_out \"/data/output_dir/$(basename "${TRIMMED_MERGED}")\""
fi
if [ -n "${FAILED_READS}" ]; then
  FASTP_COMMAND_STRING+=" --failed_out \"/data/output_dir/$(basename "${FAILED_READS}")\""
fi
FASTP_COMMAND_STRING+=" --json \"/data/output_dir/$(basename "${FASTP_JSON_REPORT}")\""
FASTP_COMMAND_STRING+=" --html \"/data/output_dir/$(basename "${FASTP_HTML_REPORT}")\""
FASTP_COMMAND_STRING+=" --thread ${NUM_THREADS}"
FASTP_COMMAND_STRING+=" --detect_adapter_for_pe"
FASTP_COMMAND_STRING+=" --cut_front"
FASTP_COMMAND_STRING+=" --cut_tail"
FASTP_COMMAND_STRING+=" --cut_tail_window_size 10"
FASTP_COMMAND_STRING+=" --cut_tail_mean_quality 20"
FASTP_COMMAND_STRING+=" --length_required 100"
FASTP_COMMAND_STRING+=" --trim_poly_g"
FASTP_COMMAND_STRING+=" --trim_poly_x"

# Execute fastp inside the Apptainer container using bash -c for Conda activation.
# Bind mounts ensure data is accessible within the container.
# We bind individual input files for clarity and explicit access.
apptainer exec \
    --bind "${INPUT_DIR}":/data/input_dir:ro \
    --bind "${OUTPUT_SAMPLE_DIR}":/data/output_dir \
    "${FASTP_CONTAINER}" \
    bash -c " \
        source \"${CONDA_HOME}/etc/profile.d/conda.sh\" && \
        conda activate fastp_env && \
        ${FASTP_COMMAND_STRING} \
    "
FASTP_EXIT_CODE=$?

#------------------------------------------------------------------------------#
#                           5. Error Handling & Completion                     #
#------------------------------------------------------------------------------#
if [ ${FASTP_EXIT_CODE} -ne 0 ]; then
    echo "ERROR: fastp execution failed with exit code ${FASTP_EXIT_CODE}."
    exit ${FASTP_EXIT_CODE}
fi

echo "----------------------------------------------------------------"
echo "Fastp QC/Trimming completed successfully for ${METAGENOME_ID}."
echo "Output trimmed reads: ${TRIMMED_R1}, ${TRIMMED_R2}"
echo "Reports: ${FASTP_JSON_REPORT}, ${FASTP_HTML_REPORT}"
echo "Job finished at: $(date)"
echo "----------------------------------------------------------------"

exit 0
