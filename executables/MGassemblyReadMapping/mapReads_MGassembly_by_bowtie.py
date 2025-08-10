#!/usr/bin/env python3
import argparse
import subprocess
import os
import sys
import logging

# --- Setup Logging ---
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s', stream=sys.stderr)
logger = logging.getLogger(__name__)

def run_mapping(
    bowtie2_index_folder: str = None,
    r1_reads: list = None,
    r2_reads: list = None,
    single_reads: list = None,
    output_dir: str = None,
    assembly_name: str = None,
    metagenome_name: str = None,
    threads: int = None
):
    """
    Constructs and executes the bowtie2 command.

    Args:
        bowtie2_index_folder (str): Directory for bowtie2 index files. (Required)
        r1_reads (list): List of paths to R1 paired-end FASTQ files. (Required)
        r2_reads (list): List of paths to R2 paired-end FASTQ files. (Required)
        single_reads (list, optional): List of paths to single-end FASTQ files.
                                       These will be passed with -s.
        output_dir (str): Base directory for mapping output. (Required)
        assembly_name (str): A name for this specific assembly, used to create
                             a subdirectory within output_dir. (Required)
        metagenome_name (str): A name for the metagenome being mapped to the assembly,
                               used to create the output file name. (Required)
        threads (int): Number of threads to use for bowtie. (Required)
    """

    # --- Input Validation (existing) ---
    if not bowtie2_index_folder:
        logger.error("Error: Index location is mandatory.")
        sys.exit(1)
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
    if not metagenome_name:
        logger.error("Error: Metagenome name is mandatory.")
        sys.exit(1)
    if threads is None:
        logger.error("Error: Number of threads is mandatory.")
        sys.exit(1)

    # Construct the full output path
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Mapping output will be written to: {output_dir}")

    # Construct the full index prefix path
    index_prefix = os.path.join(bowtie2_index_folder, f"{assembly_name}_bowtie2_index")

    # Start building the command
    cmd = ["bowtie2"]
    cmd.extend(["-x", index_prefix])

    # Add forward, reverse, and single reads to command
    r1_files_str = ",".join(r1_reads)
    r2_files_str = ",".join(r2_reads)
    single_reads_str = ",".join(single_reads) if single_reads else ""
    cmd.extend(["-1", r1_files_str,
                "-2", r2_files_str])
    if single_reads_str:
        cmd.extend(["-U", single_reads_str])
    
    # Add thread number
    cmd.extend(["-p", str(threads)])

    # Add output location
    output_file = os.path.join(output_dir, f"{metagenome_name}_to_{assembly_name}_bowtie2.sam")
    cmd.extend(["-S", output_file])

    logger.info(f"Executing Bowtie2 mapping command: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True, text=True, capture_output=False)
        logger.info("Bowtie2 execution completed successfully.")
    except subprocess.CalledProcessError as e:
        logger.error(f"Bowtie2 command failed with exit code {e.returncode}")
        logger.error(f"Command: {e.cmd}")
        sys.exit(e.returncode)
    except FileNotFoundError:
        logger.error("Error: 'bowtie2' command not found. "
                     "Ensure Bowtie2 is installed and in your PATH within the container.")
        sys.exit(1)
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}")
        sys.exit(1)

def samtools_view(
    output_dir: str = None,
    assembly_name: str = None,
    metagenome_name: str = None
    ):
    """
    Converts the SAM file to BAM format using samtools.

    Args:
        output_dir (str): Base directory for mapping output. (Required)
        assembly_name (str): A name for this specific assembly, used in the naming
                             of the mapping file. (Required)
        metagenome_name (str): A name for the metagenome being mapped to the assembly,
                               used to create the output file name. (Required)
    """

    # --- Input Validation (existing) ---
    if not output_dir:
        logger.error("Error: Output directory is mandatory for samtools conversion.")
        sys.exit(1)
    if not assembly_name:
        logger.error("Error: Assembly name is mandatory for samtools conversion.")
        sys.exit(1)
    if not metagenome_name:
        logger.error("Error: Metagenome name is mandatory for samtools conversion.")
        sys.exit(1)
    
    # Construct the full input SAM file path
    sam_file = os.path.join(output_dir, f"{metagenome_name}_to_{assembly_name}_bowtie2.sam")
    # Construct the full output BAM file path
    bam_file = os.path.join(output_dir, f"{metagenome_name}_to_{assembly_name}_bowtie2_unsorted.bam")
    
    # Set up the samtools command
    cmd = ["samtools", "view", "-b", sam_file, "-o", bam_file]
    
    logger.info(f"Converting SAM file {sam_file} to BAM format at {bam_file}: {' '.join(cmd)}")

def samtools_sort(
    output_dir: str = None,
    assembly_name: str = None,
    metagenome_name: str = None
    ):
    """
    Sorts the BAM file using samtools.

    Args:
        output_dir (str): Base directory for mapping output. (Required)
        assembly_name (str): A name for this specific assembly, used in the naming
                             of the mapping file. (Required)
        metagenome_name (str): A name for the metagenome being mapped to the assembly,
                               used to create the output file name. (Required)
    """

    # --- Input Validation (existing) ---
    if not output_dir:
        logger.error("Error: Output directory is mandatory for samtools sorting.")
        sys.exit(1)
    if not assembly_name:
        logger.error("Error: Assembly name is mandatory for samtools sorting.")
        sys.exit(1)
    if not metagenome_name:
        logger.error("Error: Metagenome name is mandatory for samtools sorting.")
        sys.exit(1)

    # Construct the full input BAM file path
    bam_file = os.path.join(output_dir, f"{metagenome_name}_to_{assembly_name}_bowtie2_unsorted.bam")
    # Construct the full output sorted BAM file path
    sorted_bam_file = os.path.join(output_dir, f"{metagenome_name}_to_{assembly_name}_bowtie2.bam")

    # Set up the samtools sort command
    cmd = ["samtools", "sort"]
    cmd.extend(["-@", str(threads)])
    cmd.extend([bam_file, "-o", sorted_bam_file])

    logger.info(f"Sorting BAM file {bam_file} to {sorted_bam_file}: {' '.join(cmd)}")

def samtools_index(
    output_dir: str = None,
    assembly_name: str = None,
    metagenome_name: str = None
    ):
    """
    Indexes the sorted BAM file using samtools.

    Args:
        output_dir (str): Base directory for mapping output. (Required)
        assembly_name (str): A name for this specific assembly, used in the naming
                             of the mapping file. (Required)
        metagenome_name (str): A name for the metagenome being mapped to the assembly,
                               used to create the output file name. (Required)
    """

    # --- Input Validation (existing) ---
    if not output_dir:
        logger.error("Error: Output directory is mandatory for samtools indexing.")
        sys.exit(1)
    if not assembly_name:
        logger.error("Error: Assembly name is mandatory for samtools indexing.")
        sys.exit(1)
    if not metagenome_name:
        logger.error("Error: Metagenome name is mandatory for samtools indexing.")
        sys.exit(1)

    # Construct the full input sorted BAM file path
    sorted_bam_file = os.path.join(output_dir, f"{metagenome_name}_to_{assembly_name}_bowtie2.bam")

    # Set up the samtools index command
    cmd = ["samtools", "index", sorted_bam_file]

    logger.info(f"Indexing BAM file {sorted_bam_file}: {' '.join(cmd)}")


def main():
    parser = argparse.ArgumentParser(
        description="Wrapper script for Bowtie mapping within Apptainer containers.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "--bowtie2_index_folder",
        required=True,
        help="Directory where bowtie2 index files are located. Required.\n"
             "Example: --bowtie2_index_folder /path/to/index_directory"
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
        nargs='+',
        default=[],
        help="Path(s) to single-end FASTQ.GZ file(s). Optional.\n"
             "If multiple files are to be provided, separate them with spaces.\n"
             "These are typically fastp's --pe-s outputs.\n"
             "Example: --single R1_unpaired.fastq.gz R2_unpaired.fastq.gz"
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
        help="Name of the assembly being mapped to the metagenome.\n"
                "This will be used to create a subdirectory within the output directory.\n"
                "Also used to create the output file name for the mapping results."
    )
    parser.add_argument(
        "--metagenome-name",
        type=str,
        required=True,
        help="Name of the metagenome being mapped to the assembly.\n"
             "This will be used to create the output file name."
    )
    parser.add_argument(
        "--threads",
        type=int,
        required=True,
        help="Number of threads (CPUs) to use for bowtie2."
    )

    args = parser.parse_args()

    # Pass the arguments to the main function
    run_mapping(
        bowtie2_index_folder=args.bowtie2_index_folder,
        r1_reads=args.r1,
        r2_reads=args.r2,
        single_reads=args.single,
        output_dir=args.output_dir,
        assembly_name=args.assembly_name,
        metagenome_name=args.metagenome_name,
        threads=args.threads
    )
    # Call samtools functions
    samtools_view(
        output_dir=args.output_dir,
        assembly_name=args.assembly_name,
        metagenome_name=args.metagenome_name
    )
    samtools_sort(
        output_dir=args.output_dir,
        assembly_name=args.assembly_name,
        metagenome_name=args.metagenome_name
    )
    samtools_index(
        output_dir=args.output_dir,
        assembly_name=args.assembly_name,
        metagenome_name=args.metagenome_name
    )


if __name__ == "__main__":
    main()