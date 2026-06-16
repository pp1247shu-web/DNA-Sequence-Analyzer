"""
main.py
-------
DNA Sequence Analyzer — Entry Point
====================================
Interactive command-line tool for bioinformatics sequence analysis.

Author  : [Your Name]
Version : 1.0.0
Python  : 3.10+

Usage
-----
    python main.py

The program offers a menu-driven interface to:
  • Analyse sequences typed directly by the user
  • Parse sequences from FASTA (.fasta / .fa) or text files
  • Generate charts, text reports, and CSV exports automatically
"""

import os
import sys
from pathlib import Path

# ── Make src/ importable when running from the project root ──────────────────
sys.path.insert(0, str(Path(__file__).parent / "src"))

from analyzer   import analyze_sequence, validate_sequence, sanitize_sequence
from parser     import load_sequences
from visualizer import plot_base_counts, plot_gc_comparison
from reporter   import generate_text_report, generate_csv_report


# --------------------------------------------------------------------------- #
#  Display helpers
# --------------------------------------------------------------------------- #

def print_banner() -> None:
    print("""
╔══════════════════════════════════════════════════════════════╗
║           DNA SEQUENCE ANALYZER  v1.0.0                     ║
║     Bioinformatics Portfolio Project — Python                ║
╚══════════════════════════════════════════════════════════════╝
    """)


def print_result_summary(result: dict) -> None:
    """Pretty-print a single analysis result to stdout."""
    bc   = result["base_counts"]
    freq = result["nucleotide_frequency"]
    rc   = result["reverse_complement"]
    rc_display = rc if len(rc) <= 60 else f"{rc[:30]}...{rc[-30:]}"

    print(f"""
  ┌─ {result['name']} {"─" * max(0, 50 - len(result['name']))}
  │  Length        : {result['length']:,} bp
  │  Counts        :  A={bc['A']}  T={bc['T']}  G={bc['G']}  C={bc['C']}
  │  Frequency     :  A={freq['A']}%  T={freq['T']}%  G={freq['G']}%  C={freq['C']}%
  │  GC content    : {result['gc_content']} %
  │  AT content    : {result['at_content']} %
  │  Most common   : {result['most_common']['base']} ({result['most_common']['count']:,}×)
  │  Rev. compl.   : {rc_display}
  └{"─" * 55}""")


def print_menu() -> None:
    print("""
  ┌─────────────────────────────────┐
  │         MAIN MENU               │
  │  1. Analyse manual sequence     │
  │  2. Analyse from file (FASTA/TXT│
  │  3. Load sample FASTA           │
  │  4. Exit                        │
  └─────────────────────────────────┘""")


# --------------------------------------------------------------------------- #
#  Analysis workflow
# --------------------------------------------------------------------------- #

def run_analysis(sequences: list[tuple[str, str]], source_label: str) -> None:
    """
    Given a list of (name, raw_sequence) tuples:
      1. Validate each sequence
      2. Run full analysis
      3. Print summary to terminal
      4. Save charts, text report, and CSV
    """
    results = []
    skipped = 0

    for name, raw_seq in sequences:
        try:
            sequence = sanitize_sequence(raw_seq)
        except ValueError as exc:
            print(f"\n  [SKIP] {name}: {exc}")
            skipped += 1
            continue

        is_valid, invalid_chars = validate_sequence(sequence)
        if not is_valid:
            print(
                f"\n  [WARNING] '{name}' contains invalid characters: "
                f"{invalid_chars}. Skipping."
            )
            skipped += 1
            continue

        result = analyze_sequence(name, sequence)
        results.append(result)
        print_result_summary(result)

    if not results:
        print("\n  No valid sequences to process.")
        return

    print(f"\n  ✓ Analysed {len(results)} sequence(s)  |  Skipped: {skipped}")

    # ── Charts ──────────────────────────────────────────────────────────────
    print("\n  Generating charts...")
    chart_paths = []
    for r in results:
        path = plot_base_counts(r, output_dir="charts")
        chart_paths.append(path)
        print(f"    • Bar chart saved: {path}")

    if len(results) > 1:
        gc_chart = plot_gc_comparison(results, output_dir="charts")
        print(f"    • GC comparison chart saved: {gc_chart}")

    # ── Reports ─────────────────────────────────────────────────────────────
    print("\n  Saving reports...")
    txt_path = generate_text_report(results, output_dir="reports",
                                    source_label=source_label)
    csv_path = generate_csv_report(results, output_dir="reports")

    print(f"    • Text report : {txt_path}")
    print(f"    • CSV export  : {csv_path}")
    print("\n  ✓ All done!\n")


# --------------------------------------------------------------------------- #
#  Menu handlers
# --------------------------------------------------------------------------- #

def handle_manual_input() -> None:
    print("\n  Enter your DNA sequence (A, T, G, C only).")
    print("  You may optionally provide a name for the sequence.")
    name = input("  Sequence name (press Enter to skip): ").strip() or "Manual_Input"
    seq  = input("  DNA sequence: ").strip()

    if not seq:
        print("  [ERROR] No sequence entered.")
        return

    run_analysis([(name, seq)], source_label=f"Manual input: {name}")


def handle_file_input() -> None:
    filepath = input("\n  Enter path to FASTA or TXT file: ").strip()

    try:
        sequences = load_sequences(filepath)
        print(f"\n  Loaded {len(sequences)} sequence(s) from '{filepath}'.")
        run_analysis(sequences, source_label=filepath)
    except (FileNotFoundError, ValueError) as exc:
        print(f"\n  [ERROR] {exc}")


def handle_sample_fasta() -> None:
    sample_path = Path(__file__).parent / "data" / "sample_sequences.fasta"
    if not sample_path.exists():
        print(f"\n  [ERROR] Sample file not found at: {sample_path}")
        print("  Please ensure 'data/sample_sequences.fasta' exists.")
        return

    try:
        sequences = load_sequences(str(sample_path))
        print(f"\n  Loaded {len(sequences)} sequence(s) from sample FASTA.")
        run_analysis(sequences, source_label=str(sample_path))
    except (FileNotFoundError, ValueError) as exc:
        print(f"\n  [ERROR] {exc}")


# --------------------------------------------------------------------------- #
#  Entry point
# --------------------------------------------------------------------------- #

def main() -> None:
    print_banner()

    menu_actions = {
        "1": handle_manual_input,
        "2": handle_file_input,
        "3": handle_sample_fasta,
        "4": None,    # exit
    }

    while True:
        print_menu()
        choice = input("  Select an option (1–4): ").strip()

        if choice == "4":
            print("\n  Goodbye! Good luck with your MSc applications.\n")
            sys.exit(0)

        action = menu_actions.get(choice)
        if action is None:
            print("  [ERROR] Invalid option. Please enter 1, 2, 3, or 4.")
        else:
            action()


if __name__ == "__main__":
    main()
