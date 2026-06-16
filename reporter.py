"""
reporter.py
-----------
Generates structured plain-text and (optionally) CSV reports
from DNA analysis results.  Reports are timestamped and saved
to the /reports directory.
"""

import csv
import os
from datetime import datetime
from pathlib import Path


SEPARATOR = "=" * 70
THIN_SEP  = "-" * 70


def _ensure_dir(directory: str) -> Path:
    p = Path(directory)
    p.mkdir(parents=True, exist_ok=True)
    return p


# --------------------------------------------------------------------------- #
#  Text report
# --------------------------------------------------------------------------- #

def _format_single(result: dict) -> str:
    """Format one analysis result as a readable text block."""
    bc   = result["base_counts"]
    freq = result["nucleotide_frequency"]
    mc   = result["most_common"]
    rc   = result["reverse_complement"]

    # Show only first/last 40 nt for very long sequences
    seq = result["sequence"]
    if len(seq) > 80:
        seq_display = f"{seq[:40]} ... {seq[-40:]}  [{len(seq)} nt total]"
    else:
        seq_display = seq

    rc_display = rc if len(rc) <= 80 else f"{rc[:40]} ... {rc[-40:]}"

    lines = [
        SEPARATOR,
        f"  SEQUENCE  : {result['name']}",
        THIN_SEP,
        f"  Input     : {seq_display}",
        f"  Length    : {result['length']:,} bp",
        "",
        "  NUCLEOTIDE COUNTS",
        f"    A : {bc['A']:>8,}   ({freq['A']:>7.4f} %)",
        f"    T : {bc['T']:>8,}   ({freq['T']:>7.4f} %)",
        f"    G : {bc['G']:>8,}   ({freq['G']:>7.4f} %)",
        f"    C : {bc['C']:>8,}   ({freq['C']:>7.4f} %)",
        "",
        "  CONTENT",
        f"    GC content : {result['gc_content']:>7.4f} %",
        f"    AT content : {result['at_content']:>7.4f} %",
        "",
        "  MOST COMMON NUCLEOTIDE",
        f"    {mc['base']}  →  {mc['count']:,} occurrences",
        "",
        "  REVERSE COMPLEMENT (5'→3')",
        f"    {rc_display}",
        SEPARATOR,
        "",
    ]
    return "\n".join(lines)


def generate_text_report(
    results: list[dict],
    output_dir: str = "reports",
    source_label: str = "User / File Input"
) -> str:
    """
    Write a complete analysis report to a timestamped .txt file.

    Parameters
    ----------
    results      : list of result dicts from analyzer.analyze_sequence()
    output_dir   : directory to save the report
    source_label : label for the data source shown in the header

    Returns
    -------
    str : absolute path to the saved report file
    """
    _ensure_dir(output_dir)

    timestamp   = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    file_ts     = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_path = Path(output_dir) / f"dna_report_{file_ts}.txt"

    header = [
        SEPARATOR,
        "        DNA SEQUENCE ANALYZER — ANALYSIS REPORT",
        SEPARATOR,
        f"  Generated  : {timestamp}",
        f"  Source     : {source_label}",
        f"  Sequences  : {len(results)}",
        SEPARATOR,
        "",
    ]

    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(header) + "\n")
        for result in results:
            fh.write(_format_single(result))

        fh.write(SEPARATOR + "\n")
        fh.write("  END OF REPORT\n")
        fh.write(SEPARATOR + "\n")

    return str(report_path.resolve())


# --------------------------------------------------------------------------- #
#  CSV export (machine-readable, great for downstream analysis)
# --------------------------------------------------------------------------- #

def generate_csv_report(
    results: list[dict],
    output_dir: str = "reports"
) -> str:
    """
    Export key statistics as a CSV file for further analysis in
    pandas, R, or Excel.
    """
    _ensure_dir(output_dir)

    file_ts  = datetime.now().strftime("%Y%m%d_%H%M%S")
    csv_path = Path(output_dir) / f"dna_stats_{file_ts}.csv"

    fieldnames = [
        "name", "length",
        "A", "T", "G", "C",
        "A_pct", "T_pct", "G_pct", "C_pct",
        "gc_content", "at_content",
        "most_common_base"
    ]

    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()

        for r in results:
            writer.writerow({
                "name":             r["name"],
                "length":           r["length"],
                "A":                r["base_counts"]["A"],
                "T":                r["base_counts"]["T"],
                "G":                r["base_counts"]["G"],
                "C":                r["base_counts"]["C"],
                "A_pct":            r["nucleotide_frequency"]["A"],
                "T_pct":            r["nucleotide_frequency"]["T"],
                "G_pct":            r["nucleotide_frequency"]["G"],
                "C_pct":            r["nucleotide_frequency"]["C"],
                "gc_content":       r["gc_content"],
                "at_content":       r["at_content"],
                "most_common_base": r["most_common"]["base"],
            })

    return str(csv_path.resolve())
