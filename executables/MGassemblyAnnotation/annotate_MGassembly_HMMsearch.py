#!/usr/bin/env python3
import argparse
import logging
import os
import subprocess
import sys


# --- Setup Logging ---
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s', stream=sys.stderr)
logger = logging.getLogger(__name__)


def parse_hmmer_domtblout(domtblout_path: str) -> set:
    """
    Parse HMMER --domtblout output and return a set of matched target sequence IDs.
    """
    hit_ids = set()

    with open(domtblout_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if not fields:
                continue
            # In hmmsearch domtblout, the first column is target name (sequence ID).
            hit_ids.add(fields[0])

    return hit_ids


def parse_hmmer_tblout(tblout_path: str) -> set:
    """
    Parse HMMER --tblout output and return a set of matched target sequence IDs.
    """
    hit_ids = set()

    with open(tblout_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if not fields:
                continue
            # In hmmsearch tblout, the first column is target name (sequence ID).
            hit_ids.add(fields[0])

    return hit_ids


def parse_domtblout_for_coverage(
    domtblout_path: str,
    min_hmm_coverage: float = None,
    min_target_coverage: float = None,
) -> set:
    """
    Parse HMMER --domtblout and return target IDs passing optional coverage filters.

    Coverage definitions:
    - HMM coverage: (hmm_to - hmm_from + 1) / qlen
    - Target coverage: (ali_to - ali_from + 1) / tlen
    """
    passing_ids = set()

    with open(domtblout_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 22:
                continue

            target_id = fields[0]

            try:
                tlen = float(fields[2])
                qlen = float(fields[5])
                hmm_from = float(fields[15])
                hmm_to = float(fields[16])
                ali_from = float(fields[17])
                ali_to = float(fields[18])
            except ValueError:
                continue

            hmm_cov = (hmm_to - hmm_from + 1.0) / qlen if qlen > 0 else 0.0
            target_cov = (ali_to - ali_from + 1.0) / tlen if tlen > 0 else 0.0

            if min_hmm_coverage is not None and hmm_cov < min_hmm_coverage:
                continue
            if min_target_coverage is not None and target_cov < min_target_coverage:
                continue

            passing_ids.add(target_id)

    return passing_ids


def write_hit_fasta(input_faa: str, output_faa: str, hit_ids: set) -> int:
    """
    Write FASTA entries from input_faa whose sequence IDs are in hit_ids.
    Returns the number of sequences written.
    """
    if not hit_ids:
        # Still create an empty output file for consistent downstream behavior.
        open(output_faa, "w", encoding="utf-8").close()
        return 0

    written = 0
    current_header = None
    current_id = None
    current_seq_lines = []

    def flush_record(out_handle):
        nonlocal written
        if current_header is None:
            return
        if current_id in hit_ids:
            out_handle.write(current_header + "\n")
            for seq_line in current_seq_lines:
                out_handle.write(seq_line)
            if not current_seq_lines or not current_seq_lines[-1].endswith("\n"):
                out_handle.write("\n")
            written += 1

    with open(input_faa, "r", encoding="utf-8") as in_handle, open(output_faa, "w", encoding="utf-8") as out_handle:
        for line in in_handle:
            if line.startswith(">"):
                flush_record(out_handle)
                current_header = line.rstrip("\n")
                current_id = current_header[1:].split()[0] if current_header[1:].strip() else ""
                current_seq_lines = []
            else:
                current_seq_lines.append(line)

        flush_record(out_handle)

    return written


def run_hmmsearch(
    hmm_model: str,
    input_orf_files: list,
    output_dir: str,
    evalue_cutoff: float = None,
    cut_tc: bool = False,
    cut_ga: bool = False,
    cut_nc: bool = False,
    dom_evalue_cutoff: float = None,
    min_hmm_coverage: float = None,
    min_target_coverage: float = None,
    cpu: int = 1,
    hmm_stdout_ext: str = ".hmmsearch.txt",
    tblout_ext: str = ".hmmsearch.tblout",
    domtblout_ext: str = ".hmmsearch.domtblout",
    hits_faa_ext: str = ".hmm_hits.faa",
    output_basename: str = None,
    keep_domtblout: bool = False,
):
    """
    Run HMMER hmmsearch in a container for each ORF amino acid file and extract matching sequences.
    """
    if not os.path.exists(hmm_model):
        logger.error(f"Error: HMM model file does not exist: {hmm_model}")
        sys.exit(1)
    if not input_orf_files:
        logger.error("Error: At least one ORF amino acid input file is required.")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Output directory: {output_dir}")

    for input_faa in input_orf_files:
        if not os.path.exists(input_faa):
            logger.error(f"Error: Input ORF amino acid file does not exist: {input_faa}")
            sys.exit(1)

        orf_base = os.path.splitext(os.path.basename(input_faa))[0]
        hmm_base = os.path.splitext(os.path.basename(hmm_model))[0]
        if output_basename:
            base_name = output_basename if len(input_orf_files) == 1 else f"{output_basename}__{orf_base}"
        else:
            base_name = f"{orf_base}__{hmm_base}"

        hmm_stdout_path = os.path.join(output_dir, f"{base_name}{hmm_stdout_ext}")
        tblout_path = os.path.join(output_dir, f"{base_name}{tblout_ext}")
        domtblout_path = os.path.join(output_dir, f"{base_name}{domtblout_ext}")
        hits_faa_path = os.path.join(output_dir, f"{base_name}{hits_faa_ext}")

        needs_domtblout = keep_domtblout or min_hmm_coverage is not None or min_target_coverage is not None

        cmd = [
            "hmmsearch",
            "--cpu",
            str(cpu),
            "--tblout",
            tblout_path,
        ]

        if needs_domtblout:
            cmd.extend(["--domtblout", domtblout_path])

        if cut_tc:
            cmd.append("--cut_tc")
        elif cut_ga:
            cmd.append("--cut_ga")
        elif cut_nc:
            cmd.append("--cut_nc")
        else:
            cmd.extend(["-E", str(evalue_cutoff)])

        if dom_evalue_cutoff is not None:
            cmd.extend(["--domE", str(dom_evalue_cutoff)])

        cmd.extend([hmm_model, input_faa])

        logger.info(f"Processing ORF file: {input_faa}")
        logger.info(f"HMMER stdout output: {hmm_stdout_path}")
        logger.info(f"HMMER tblout output: {tblout_path}")
        if needs_domtblout:
            logger.info(f"HMMER domtblout output: {domtblout_path}")
        logger.info(f"Extracted hit FASTA output: {hits_faa_path}")
        logger.info(f"Executing command: {' '.join(cmd)}")

        try:
            with open(hmm_stdout_path, "w", encoding="utf-8") as out_handle:
                subprocess.run(cmd, check=True, text=True, stdout=out_handle, stderr=sys.stderr)
        except subprocess.CalledProcessError as e:
            logger.error(f"hmmsearch command failed with exit code {e.returncode}")
            logger.error(f"Command: {e.cmd}")
            sys.exit(e.returncode)
        except FileNotFoundError:
            logger.error(
                "Error: hmmsearch command not found. "
                "Ensure HMMER is installed and available in PATH (e.g., via apptainer exec)."
            )
            sys.exit(1)
        except Exception as e:
            logger.error(f"An unexpected error occurred while running hmmsearch: {e}")
            sys.exit(1)

        seq_hit_ids = parse_hmmer_tblout(tblout_path)

        if min_hmm_coverage is not None or min_target_coverage is not None:
            dom_hit_ids = parse_domtblout_for_coverage(
                domtblout_path=domtblout_path,
                min_hmm_coverage=min_hmm_coverage,
                min_target_coverage=min_target_coverage,
            )
            hit_ids = seq_hit_ids.intersection(dom_hit_ids)
            logger.info(f"Sequence-level hits from tblout: {len(seq_hit_ids)}")
            logger.info(f"Hits passing coverage filters: {len(hit_ids)}")
        else:
            hit_ids = seq_hit_ids

        n_written = write_hit_fasta(input_faa=input_faa, output_faa=hits_faa_path, hit_ids=hit_ids)
        logger.info(f"Matched sequence IDs retained: {len(hit_ids)}")
        logger.info(f"Amino acid sequences written to hit FASTA: {n_written}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Wrapper script for HMM-based ORF annotation using hmmsearch (expects hmmsearch in PATH).",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "--hmm-model",
        type=str,
        required=True,
        help="Path to the HMM profile/model file for the protein of interest. Required.",
    )
    parser.add_argument(
        "--input-aa-orfs",
        nargs="+",
        required=True,
        help="Path(s) to amino acid ORF FASTA file(s) to search. Required.",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Output directory for hmmsearch results and extracted hit sequences. Required.",
    )
    parser.add_argument(
        "--cpu",
        type=int,
        default=1,
        help="Number of CPU threads for hmmsearch (--cpu). Default: 1",
    )

    cutoff_group = parser.add_mutually_exclusive_group(required=True)
    cutoff_group.add_argument(
        "--evalue-cutoff",
        type=float,
        help="Sequence-level E-value cutoff passed to hmmsearch as -E.",
    )
    cutoff_group.add_argument(
        "--cut-tc",
        action="store_true",
        help="Use trusted cutoff from model annotation (hmmsearch --cut_tc).",
    )
    cutoff_group.add_argument(
        "--cut-ga",
        action="store_true",
        help="Use curated gathering cutoff from model annotation (hmmsearch --cut_ga).",
    )
    cutoff_group.add_argument(
        "--cut-nc",
        action="store_true",
        help="Use curated noise cutoff from model annotation (hmmsearch --cut_nc).",
    )

    parser.add_argument(
        "--dom-evalue-cutoff",
        type=float,
        default=None,
        help="Optional domain-level E-value cutoff passed as --domE.",
    )
    parser.add_argument(
        "--min-hmm-coverage",
        type=float,
        default=None,
        help="Optional minimum HMM coverage fraction (0-1) for domain matches.",
    )
    parser.add_argument(
        "--min-target-coverage",
        type=float,
        default=None,
        help="Optional minimum target sequence coverage fraction (0-1) for domain matches.",
    )
    parser.add_argument(
        "--output-basename",
        type=str,
        default=None,
        help=(
            "Optional basename for output files. "
            "If omitted, defaults to <orf_basename>__<hmm_basename>. "
            "If multiple ORF files are supplied, each output becomes <output_basename>__<orf_basename>."
        ),
    )
    parser.add_argument(
        "--keep-domtblout",
        action="store_true",
        help="Write domtblout output even when coverage filters are not requested.",
    )
    parser.add_argument(
        "--hmm-stdout-ext",
        type=str,
        default=".hmmsearch.txt",
        help="Extension for saved hmmsearch stdout files. Default: .hmmsearch.txt",
    )
    parser.add_argument(
        "--tblout-ext",
        type=str,
        default=".hmmsearch.tblout",
        help="Extension for hmmsearch --tblout files. Default: .hmmsearch.tblout",
    )
    parser.add_argument(
        "--domtblout-ext",
        type=str,
        default=".hmmsearch.domtblout",
        help="Extension for hmmsearch --domtblout files. Default: .hmmsearch.domtblout",
    )
    parser.add_argument(
        "--hits-faa-ext",
        type=str,
        default=".hmm_hits.faa",
        help="Extension for extracted hit sequence FASTA files. Default: .hmm_hits.faa",
    )

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    if args.cpu < 1:
        parser.error("--cpu must be >= 1")
    if args.evalue_cutoff is not None and args.evalue_cutoff <= 0:
        parser.error("--evalue-cutoff must be > 0")
    if args.dom_evalue_cutoff is not None and args.dom_evalue_cutoff <= 0:
        parser.error("--dom-evalue-cutoff must be > 0")
    if args.min_hmm_coverage is not None and not (0 < args.min_hmm_coverage <= 1):
        parser.error("--min-hmm-coverage must be in (0, 1]")
    if args.min_target_coverage is not None and not (0 < args.min_target_coverage <= 1):
        parser.error("--min-target-coverage must be in (0, 1]")

    run_hmmsearch(
        hmm_model=args.hmm_model,
        input_orf_files=args.input_aa_orfs,
        output_dir=args.output_dir,
        evalue_cutoff=args.evalue_cutoff,
        cut_tc=args.cut_tc,
        cut_ga=args.cut_ga,
        cut_nc=args.cut_nc,
        dom_evalue_cutoff=args.dom_evalue_cutoff,
        min_hmm_coverage=args.min_hmm_coverage,
        min_target_coverage=args.min_target_coverage,
        cpu=args.cpu,
        hmm_stdout_ext=args.hmm_stdout_ext,
        tblout_ext=args.tblout_ext,
        domtblout_ext=args.domtblout_ext,
        hits_faa_ext=args.hits_faa_ext,
        output_basename=args.output_basename,
        keep_domtblout=args.keep_domtblout,
    )


if __name__ == "__main__":
    main()