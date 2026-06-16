"""
parser.py
---------
Handles reading DNA sequences from FASTA (.fasta / .fa) files
and plain text files.

FASTA format reminder
---------------------
>SequenceID optional description
ATGCATGCATGC...
ATGCATGCATGC...   <- sequence may span multiple lines
"""

import os
from pathlib import Path


# --------------------------------------------------------------------------- #
#  FASTA parser
# --------------------------------------------------------------------------- #

def parse_fasta(filepath: str) -> list[tuple[str, str]]:
    """
    Parse a FASTA file and return a list of (header, sequence) tuples.

    Parameters
    ----------
    filepath : str
        Absolute or relative path to a .fasta / .fa file.

    Returns
    -------
    List of (identifier, sequence) tuples.

    Raises
    ------
    FileNotFoundError  – path does not exist
    ValueError         – file is empty or contains no valid FASTA records
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {filepath}")
    if path.stat().st_size == 0:
        raise ValueError(f"File is empty: {filepath}")

    records = []
    current_header = None
    current_seq_parts: list[str] = []

    with open(path, "r", encoding="utf-8") as fh:
        for line_num, line in enumerate(fh, start=1):
            line = line.strip()

            if not line or line.startswith(";"):   # skip blank lines & comments
                continue

            if line.startswith(">"):
                # Save the previous record before starting a new one
                if current_header is not None:
                    seq = "".join(current_seq_parts).upper()
                    records.append((current_header, seq))

                current_header = line[1:].strip()  # strip the ">"
                current_seq_parts = []

            else:
                if current_header is None:
                    raise ValueError(
                        f"Line {line_num}: sequence data found before any FASTA header."
                    )
                current_seq_parts.append(line)

    # Don't forget the last record
    if current_header is not None and current_seq_parts:
        seq = "".join(current_seq_parts).upper()
        records.append((current_header, seq))

    if not records:
        raise ValueError(f"No valid FASTA records found in: {filepath}")

    return records


# --------------------------------------------------------------------------- #
#  Plain text parser
# --------------------------------------------------------------------------- #

def parse_text_file(filepath: str) -> list[tuple[str, str]]:
    """
    Read a plain text file where each non-empty line is a separate sequence.
    Lines beginning with '#' are treated as comments.

    Returns a list of ("Sequence_N", sequence) tuples.
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    records = []
    with open(path, "r", encoding="utf-8") as fh:
        for idx, line in enumerate(fh, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            records.append((f"Sequence_{idx}", line.upper()))

    if not records:
        raise ValueError(f"No sequences found in: {filepath}")

    return records


# --------------------------------------------------------------------------- #
#  Auto-dispatch
# --------------------------------------------------------------------------- #

def load_sequences(filepath: str) -> list[tuple[str, str]]:
    """
    Automatically choose the right parser based on file extension.

    Supported: .fasta, .fa, .fna, .txt
    """
    ext = Path(filepath).suffix.lower()

    if ext in {".fasta", ".fa", ".fna"}:
        return parse_fasta(filepath)
    elif ext == ".txt":
        return parse_text_file(filepath)
    else:
        raise ValueError(
            f"Unsupported file type '{ext}'. "
            "Supported formats: .fasta, .fa, .fna, .txt"
        )
