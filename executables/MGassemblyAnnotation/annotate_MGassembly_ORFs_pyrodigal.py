#!/usr/bin/env python3
import argparse
import subprocess
import os
import sys
import logging

# --- Setup Logging ---
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s', stream=sys.stderr)
logger = logging.getLogger(__name__)


def run_pyrodigal(
    input_assemblies: list,
    output_dir: str = None,
    output_format: str = "gff",
    procedure: str = "anon",
    protein_ext: str = ".faa",
    gff_ext: str = ".gff"
):
    """
    Constructs and executes the pyrodigal command for ORF prediction.

    Args:
        input_assemblies (list): List of paths to assembly files (.fna files). (Required)
        output_dir (str): Base directory for pyrodigal output. (Required)
        output_format (str): Output format for annotations (gff, prodigal, csv). Default: "gff"
        procedure (str): Procedure mode (single, meta, anon). Default: "anon"
        protein_ext (str): Extension for protein output files. Default: ".faa"
        gff_ext (str): Extension for GFF output files. Default: ".gff"
    """

    # --- Input Validation ---
    if not input_assemblies:
        logger.error("Error: Input assembly files are mandatory.")
        sys.exit(1)
    if not output_dir:
        logger.error("Error: Output directory is mandatory.")
        sys.exit(1)

    # Construct the output path
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Pyrodigal output will be written to: {output_dir}")

    # Process each input assembly file
    for input_assembly in input_assemblies:
        if not os.path.exists(input_assembly):
            logger.error(f"Error: Input assembly file '{input_assembly}' does not exist.")
            sys.exit(1)

        # Derive output file names from input filename
        base_name = os.path.splitext(os.path.basename(input_assembly))[0]
        protein_output = os.path.join(output_dir, f"{base_name}{protein_ext}")
        gff_output = os.path.join(output_dir, f"{base_name}{gff_ext}")

        # Start building the command
        cmd = ["pyrodigal"]
        cmd.extend(["-i", input_assembly])
        cmd.extend(["-a", protein_output])
        cmd.extend(["-o", gff_output])
        cmd.extend(["-f", output_format])
        cmd.extend(["-p", procedure])

        logger.info(f"Processing assembly: {input_assembly}")
        logger.info(f"Protein output: {protein_output}")
        logger.info(f"GFF output: {gff_output}")
        logger.info(f"Executing pyrodigal command: {' '.join(cmd)}")

        try:
            # Use subprocess.run for cleaner handling of process execution
            subprocess.run(cmd, check=True, text=True, capture_output=False)
            logger.info(f"Pyrodigal execution completed successfully for {input_assembly}.")
        except subprocess.CalledProcessError as e:
            logger.error(f"Pyrodigal command failed with exit code {e.returncode}")
            logger.error(f"Command: {e.cmd}")
            sys.exit(e.returncode)
        except FileNotFoundError:
            logger.error("Error: 'pyrodigal' command not found. "
                         "Ensure pyrodigal is installed and in your PATH within the container.")
            sys.exit(1)
        except Exception as e:
            logger.error(f"An unexpected error occurred: {e}")
            sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Wrapper script for ORF prediction using Pyrodigal within Apptainer containers.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "--input-assemblies",
        nargs='+',  # Allows multiple assembly files
        required=True,
        help="Path(s) to assembly file(s) (.fna) to predict ORFs from. Required.\n"
             "Example: --input-assemblies assembly1.fna assembly2.fna"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Base output directory for pyrodigal results.\n"
             "Required.\n"
             "Example: --output-dir ./annotation_output"
    )
    parser.add_argument(
        "--output-format",
        type=str,
        default="gff",
        choices=["gff", "prodigal", "csv"],
        help="Output format for annotations (gff, prodigal, csv).\n"
             "Default: gff\n"
             "Example: --output-format gff"
    )
    parser.add_argument(
        "--procedure",
        type=str,
        default="anon",
        choices=["single", "meta", "anon"],
        help="Procedure mode for ORF prediction (single, meta, anon).\n"
             "Default: anon (anonymous mode for metagenomic data)\n"
             "Example: --procedure anon"
    )
    parser.add_argument(
        "--protein-ext",
        type=str,
        default=".faa",
        help="Extension for protein output files.\n"
             "Default: .faa\n"
             "Example: --protein-ext .faa"
    )
    parser.add_argument(
        "--gff-ext",
        type=str,
        default=".gff",
        help="Extension for GFF annotation output files.\n"
             "Default: .gff\n"
             "Example: --gff-ext .gff"
    )

    args = parser.parse_args()

    # Pass the arguments to the main function
    run_pyrodigal(
        input_assemblies=args.input_assemblies,
        output_dir=args.output_dir,
        output_format=args.output_format,
        procedure=args.procedure,
        protein_ext=args.protein_ext,
        gff_ext=args.gff_ext
    )


if __name__ == "__main__":
    main()