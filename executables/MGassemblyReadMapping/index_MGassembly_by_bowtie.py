#!/usr/bin/env python3
import argparse
import subprocess
import os
import sys
import logging

# --- Setup Logging ---
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s', stream=sys.stderr)
logger = logging.getLogger(__name__)



def run_bowtie_indexing(
    input_assembly: str = None,
    output_index_name: str = None,
    output_dir: str = None
):
    """
    Constructs and executes the bowtie2-build command.

    Args:
        input_assembly (str): Path to input assembly
        output_index_name (str): Name to be used for output file. Usually the assembly name.
        output_index_dir (str): Base directory and prefix for indexing output. (Required)
    """

    # --- Input Validation (existing) ---
    if not input_assembly:
        logger.error("Error: Assembly input is mandatory.")
        sys.exit(1)    
    if not output_index_name:
        logger.error("Error: Output index name is mandatory.")
        sys.exit(1)  
    if not output_dir:
        logger.error("Error: Output location is mandatory.")
        sys.exit(1)

    # Construct the output path
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Index output will be written to: {output_dir}")

    # Add output file name to output path
    output_path = output_dir + output_index_name + "_bowtie2_index"

    # Start building the command
    cmd = ["bowtie2-build"]
    cmd.extend([input_assembly])
    cmd.extend([output_path])

    logger.info(f"Executing indexing command: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True, text=True, capture_output=False)
        logger.info("Indexing execution completed successfully.")
    except subprocess.CalledProcessError as e:
        logger.error(f"Indexing command failed with exit code {e.returncode}")
        logger.error(f"Command: {e.cmd}")
        sys.exit(e.returncode)
    except FileNotFoundError:
        logger.error("Error: 'bowtie2-build' command not found. "
                     "Ensure bowtie2 is installed and in your PATH within the container.")
        sys.exit(1)
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="Wrapper script for indexing metagenome assembly using bowtie2 within Apptainer containers.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "--input-assembly",
        type=str,
        required=True,
        help="Path to assembly file to be indexed.\n"
             "Required.\n"
             "Example: --input-assembly assemblyID_assembly.fna"
    )
    parser.add_argument(
        "--output-name",
        type=str,
        required=True,
        help="Name to be used for output file. Usually the same as the assembly ID.\n"
             "Required.\n"
             "Example: --output-name assemblyID"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Base output directory for the assembly results.\n"
             "Required.\n"
             "Example: --output-name assemblyID"
    )

    args = parser.parse_args()

    # Pass the arguments to the main function
    run_bowtie_indexing(
        input_assembly=args.input_assembly,
        output_index_name=args.output_name,
        output_dir=args.output_dir,
    )

if __name__ == "__main__":
    main()
