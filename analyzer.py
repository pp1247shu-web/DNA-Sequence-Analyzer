"""
analyzer.py
-----------
Core bioinformatics logic for the DNA Sequence Analyzer.
Handles sequence validation, base counting, GC/AT content,
reverse complement generation, and frequency statistics.
"""

from collections import Counter
from typing import Optional


# IUPAC standard DNA nucleotides
VALID_BASES = frozenset("ATGCatgc")

COMPLEMENT_MAP = str.maketrans("ATGCatgc", "TACGtacg")


# --------------------------------------------------------------------------- #
#  Validation
# --------------------------------------------------------------------------- #

def validate_sequence(sequence: str) -> tuple[bool, set]:
    """
    Validate that all characters in the sequence are standard DNA bases.

    Parameters
    ----------
    sequence : str
        Raw DNA string (case-insensitive).

    Returns
    -------
    is_valid : bool
    invalid_chars : set   – empty if the sequence is clean
    """
    sequence = sequence.strip().upper()
    chars = set(sequence)
    invalid = chars - set("ATGC")
    return (len(invalid) == 0), invalid


def sanitize_sequence(sequence: str) -> str:
    """
    Strip whitespace/newlines and convert to uppercase.
    Raises ValueError if the sequence is empty after cleaning.
    """
    clean = "".join(sequence.split()).upper()
    if not clean:
        raise ValueError("Sequence is empty after removing whitespace.")
    return clean


# --------------------------------------------------------------------------- #
#  Core statistics
# --------------------------------------------------------------------------- #

def count_bases(sequence: str) -> dict[str, int]:
    """Return a dict with counts for A, T, G, C (and N for unknowns)."""
    counts = Counter(sequence)
    return {
        "A": counts.get("A", 0),
        "T": counts.get("T", 0),
        "G": counts.get("G", 0),
        "C": counts.get("C", 0),
    }


def calculate_gc_content(sequence: str) -> float:
    """
    GC content = (G + C) / total_length × 100.
    A GC content between 40–60 % is typical for most organisms.
    """
    if not sequence:
        return 0.0
    counts = count_bases(sequence)
    gc = counts["G"] + counts["C"]
    return round((gc / len(sequence)) * 100, 4)


def calculate_at_content(sequence: str) -> float:
    """AT content = 100 - GC content."""
    return round(100.0 - calculate_gc_content(sequence), 4)


def get_sequence_length(sequence: str) -> int:
    """Return the number of nucleotides in the sequence."""
    return len(sequence)


def most_common_nucleotide(sequence: str) -> tuple[str, int]:
    """Return (nucleotide, count) for the most frequent base."""
    counts = count_bases(sequence)
    base = max(counts, key=counts.get)
    return base, counts[base]


def nucleotide_frequency(sequence: str) -> dict[str, float]:
    """
    Return percentage frequency for each base.
    Useful for quick comparison between sequences of different lengths.
    """
    length = len(sequence)
    counts = count_bases(sequence)
    return {base: round((n / length) * 100, 4) for base, n in counts.items()}


# --------------------------------------------------------------------------- #
#  Sequence operations
# --------------------------------------------------------------------------- #

def reverse_complement(sequence: str) -> str:
    """
    Generate the reverse complement of a DNA sequence.

    Biology note
    ------------
    DNA is antiparallel and double-stranded.  The complement of a given
    5'→3' strand read 5'→3' is the reverse complement.
    Example: 5'-ATCG-3'  →  reverse complement: 5'-CGAT-3'
    """
    return sequence.translate(COMPLEMENT_MAP)[::-1]


def find_motif(sequence: str, motif: str) -> list[int]:
    """
    Return all 1-based start positions of a motif within the sequence.
    Useful for identifying restriction sites or primer binding regions.
    """
    motif = motif.upper()
    positions = []
    start = 0
    while True:
        pos = sequence.find(motif, start)
        if pos == -1:
            break
        positions.append(pos + 1)   # convert to 1-based biological indexing
        start = pos + 1
    return positions


# --------------------------------------------------------------------------- #
#  Full analysis pipeline
# --------------------------------------------------------------------------- #

def analyze_sequence(name: str, sequence: str) -> dict:
    """
    Run all analyses on a single sequence and return a unified result dict.

    Parameters
    ----------
    name     : sequence identifier (e.g. FASTA header)
    sequence : validated, uppercase DNA string

    Returns
    -------
    dict with all computed fields
    """
    bases   = count_bases(sequence)
    freq    = nucleotide_frequency(sequence)
    top_nt  = most_common_nucleotide(sequence)

    return {
        "name":               name,
        "sequence":           sequence,
        "length":             get_sequence_length(sequence),
        "base_counts":        bases,
        "gc_content":         calculate_gc_content(sequence),
        "at_content":         calculate_at_content(sequence),
        "nucleotide_frequency": freq,
        "most_common":        {"base": top_nt[0], "count": top_nt[1]},
        "reverse_complement": reverse_complement(sequence),
    }
