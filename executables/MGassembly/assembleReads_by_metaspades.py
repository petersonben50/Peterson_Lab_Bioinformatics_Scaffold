#!/usr/bin/env python3
import argparse
import subprocess
import os
import sys
import logging

# --- Setup Logging ---
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s', stream=sys.stderr)
logger = logging.getLogger(__name__)

def validate_kmers(kmer_str: str):
    """
    Validates the k-mer string format and individual k-mer values.
    K-mers must be odd integers less than 128.
    """
    if not kmer_str:
        return "" # Allow empty string if default is handled elsewhere or not required

    k_list = kmer_str.split(',')
    validated_kmers = []
    for k in k_list:
        try:
            k_int = int(k.strip())
            if not (k_int <= 150):
                raise ValueError(f"K-mer value '{k_int}' is out of the allowed range (under 128).")
            if k_int % 2 == 0:
                raise ValueError(f"K-mer value '{k_int}' must be an odd number.")
            validated_kmers.append(str(k_int))
        except ValueError as e:
            logger.error(f"Invalid k-mer value '{k}': {e}")
            sys.exit(1)
    
    return ",".join(validated_kmers)

def run_metaspades(
    r1_reads: list,
    r2_reads: list,
    single_reads: list = None,
    merged_reads: list = None,
    output_dir: str = None,
    assembly_name: str = None,
    threads: int = None,
    memory_limit_gb: int = None,
    kmers: str = None
):
    """
    Constructs and executes the spades.py command with --meta flag.

    Args:
        r1_reads (list): List of paths to R1 paired-end FASTQ files. (Required)
        r2_reads (list): List of paths to R2 paired-end FASTQ files. (Required)
        single_reads (list, optional): List of paths to single-end FASTQ files.
                                       These will be passed with -s.
        merged_reads (list, optional): List of paths to merged FASTQ files.
                                       These will be passed with --merged.
        output_dir (str): Base directory for SPAdes output. (Required)
        assembly_name (str): A name for this specific assembly, used to create
                             a subdirectory within output_dir. (Required)
        threads (int): Number of threads to use for SPAdes. (Required)
        memory_limit_gb (int, optional): Maximum RAM in Gigabytes for SPAdes to use (-m flag).
        kmers (str, optional): Comma-separated string of k-mer values.
                               e.g., "21,33,55,77,99,127".
    """

    # --- Input Validation (existing) ---
    if not r1_reads or not r2_reads:
        logger.error("Error: Paired-end R1 and R2 reads are mandatory.")
        sys.exit(1)
    if len(r1_reads) != len(r2_reads):
        logger.error("Error: Number of R1 files must match the number of R2 files.")
        sys.exit(1)
    if not output_dir:
        logger.error("Error: Output directory is mandatory.")
        sys.exit(1)
    if not assembly_name:
        logger.error("Error: Assembly name is mandatory.")
        sys.exit(1)
    if threads is None:
        logger.error("Error: Number of threads is mandatory.")
        sys.exit(1)

    # --- Validate K-mers ---
    validated_kmers_str = validate_kmers(kmers)
    logger.info(f"Using k-mers: {validated_kmers_str}")

    # Construct the full output path
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"SPAdes output will be written to: {output_dir}")

    # Start building the command
    cmd = ["spades.py", "--meta"] # Always use --meta for metaSPAdes mode
    cmd.extend(["-o", output_dir])
    cmd.extend(["-t", str(threads)])

    # Add memory limit if provided
    if memory_limit_gb is not None:
        cmd.extend(["-m", str(memory_limit_gb)])
        logger.info(f"Setting SPAdes memory limit to {memory_limit_gb} GB.")
    
    # Add k-mers
    cmd.extend(["-k", validated_kmers_str])

    # Add paired reads
    for i, (r1, r2) in enumerate(zip(r1_reads, r2_reads)):
        if i == 0: # First pair uses -1 and -2
            cmd.extend(["-1", r1, "-2", r2])
        else: # Subsequent pairs use --peK-1 and --peK-2
            cmd.extend([f"--pe{i}-1", r1, f"--pe{i}-2", r2])
    
    # Add single reads
    if single_reads:
        for s_read in single_reads:
            cmd.extend(["-s", s_read])

    # Add merged reads
    if merged_reads:
        for m_read in merged_reads:
            cmd.extend(["--merged", m_read])

    logger.info(f"Executing SPAdes command: {' '.join(cmd)}")

    try:
        # Use subprocess.run for cleaner handling of process execution
        subprocess.run(cmd, check=True, text=True, capture_output=False)
        logger.info("SPAdes execution completed successfully.")
    except subprocess.CalledProcessError as e:
        logger.error(f"SPAdes command failed with exit code {e.returncode}")
        logger.error(f"Command: {e.cmd}")
        sys.exit(e.returncode)
    except FileNotFoundError:
        logger.error("Error: 'spades.py' command not found. "
                     "Ensure SPAdes is installed and in your PATH within the container.")
        sys.exit(1)
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="Wrapper script for metaSPAdes assembly within Apptainer containers.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "--r1",
        nargs='+', # Allows multiple R1 files
        required=True,
        help="Path(s) to R1 paired-end FASTQ.GZ file(s). Required.\n"
             "Example: --r1 sample_R1.fastq.gz"
    )
    parser.add_argument(
        "--r2",
        nargs='+', # Allows multiple R2 files
        required=True,
        help="Path(s) to R2 paired-end FASTQ.GZ file(s). Required.\n"
             "Example: --r2 sample_R2.fastq.gz"
    )
    parser.add_argument(
        "--single",
        nargs='+', # Allows multiple single files
        default=[],
        help="Path(s) to single-end FASTQ.GZ file(s). Optional.\n"
             "These are typically fastp's --pe-s outputs.\n"
             "Example: --single R1_unpaired.fastq.gz R2_unpaired.fastq.gz"
    )
    parser.add_argument(
        "--merged",
        nargs='+', # Allows multiple merged files
        default=[],
        help="Path(s) to merged FASTQ.GZ file(s) from fastp. Optional.\n"
             "Example: --merged merged.fastq.gz"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Base output directory for the assembly results."
    )
    parser.add_argument(
        "--assembly-name",
        type=str,
        required=True,
        help="Unique name for this assembly run. Used to create a subdirectory "
             "within the output-dir."
    )
    parser.add_argument(
        "--threads",
        type=int,
        required=True,
        help="Number of threads (CPUs) to use for SPAdes."
    )
    parser.add_argument(
        "--memory",
        type=int,
        default=None,
        help="Maximum RAM in Gigabytes for SPAdes to use (-m flag). Optional.\n"
             "It's highly recommended to set this to match your Slurm allocation."
    )
    parser.add_argument(
        "--kmers",
        type=str,
        default="21,33,55,77,99,127", # Default k-mers as specified
        help="Comma-separated string of k-mer values to use for assembly (-k flag).\n"
             "K-mer values must be odd integers between 20 and 150.\n"
             "Default: 21,33,55,77,99,127"
    )

    args = parser.parse_args()

    # Pass the arguments to the main function
    run_metaspades(
        r1_reads=args.r1,
        r2_reads=args.r2,
        single_reads=args.single,
        merged_reads=args.merged,
        output_dir=args.output_dir,
        assembly_name=args.assembly_name,
        threads=args.threads,
        memory_limit_gb=args.memory,
        kmers=args.kmers # Pass the k-mers argument
    )

if __name__ == "__main__":
    main()