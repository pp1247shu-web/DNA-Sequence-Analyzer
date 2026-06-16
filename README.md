# DNA Sequence Analyzer

> A bioinformatics command-line tool for DNA sequence analysis, GC/AT content profiling, nucleotide frequency statistics, and automated report generation.

![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![Status](https://img.shields.io/badge/Status-Portfolio%20Project-orange)

---

## Overview

DNA Sequence Analyzer is a Python-based bioinformatics tool designed to process DNA sequences from manual input or standard FASTA/text files. It computes common sequence statistics, generates publication-quality charts, and exports structured reports — covering key skills expected in MSc Bioinformatics programmes and entry-level research roles.

---

## Features

| Feature | Description |
|---|---|
| **FASTA / TXT parsing** | Multi-record FASTA support with robust error handling |
| **Base counting** | A, T, G, C absolute counts per sequence |
| **GC / AT content** | Percentage with 4 decimal precision |
| **Nucleotide frequency** | Per-base percentage frequency |
| **Reverse complement** | Biologically accurate 5'→3' reverse complement |
| **Motif search** | 1-based positional search for any k-mer |
| **Visualisations** | Bar chart (counts + pie chart) per sequence; GC comparison chart |
| **Text reports** | Timestamped, structured `.txt` reports |
| **CSV export** | Machine-readable output for downstream analysis in pandas / R |
| **Unit tests** | Full pytest suite for core logic |

---

## Project Structure

```
dna_analyzer/
├── main.py                       # Entry point (interactive CLI)
├── requirements.txt
├── README.md
│
├── src/
│   ├── analyzer.py               # Core bioinformatics logic
│   ├── parser.py                 # FASTA / TXT file parsing
│   ├── visualizer.py             # matplotlib chart generation
│   └── reporter.py               # Text & CSV report generation
│
├── data/
│   └── sample_sequences.fasta    # 5 real biological sequences
│
├── charts/                       # Auto-generated PNG charts
├── reports/                      # Auto-generated .txt and .csv reports
└── tests/
    └── test_analyzer.py          # pytest unit tests
```

---

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/dna-sequence-analyzer.git
cd dna-sequence-analyzer

# (Recommended) Create a virtual environment
python -m venv venv
source venv/bin/activate        # macOS / Linux
venv\Scripts\activate           # Windows

# Install dependencies
pip install -r requirements.txt
```

---

## Usage

```bash
python main.py
```

An interactive menu will appear:

```
  ┌─────────────────────────────────┐
  │  1. Analyse manual sequence     │
  │  2. Analyse from file           │
  │  3. Load sample FASTA           │
  │  4. Exit                        │
  └─────────────────────────────────┘
```

---

## Sample Output

```
  ┌─ NC_000913.3 Escherichia coli K-12 ─────────────────────────
  │  Length        : 300 bp
  │  Counts        :  A=62  T=58  G=91  C=89
  │  Frequency     :  A=20.67%  T=19.33%  G=30.33%  C=29.67%
  │  GC content    : 60.0 %
  │  AT content    : 40.0 %
  │  Most common   : G (91×)
  │  Rev. compl.   : CGAT...GCTA
  └─────────────────────────────────────────────────────────────
```

---

## Running Tests

```bash
python -m pytest tests/ -v
```

---

## Biological Background

- **GC content** influences DNA melting temperature (Tm), chromatin structure, and codon usage bias. Bacterial genomes range from ~25 % (Mycoplasma) to ~75 % (Streptomyces).
- **Reverse complement** is essential for primer design, restriction enzyme analysis, and understanding antisense strands.
- **Nucleotide frequency** deviations can indicate CpG islands, AT-rich regulatory regions, or sequencing artefacts.

---

## Future Upgrades (Roadmap)

| Priority | Feature |
|---|---|
| ⭐ High | Protein translation (ORF detection, codon table) |
| ⭐ High | BLAST API integration for sequence similarity search |
| ⭐ High | CpG island detection |
| Medium | Multiple sequence alignment (Biopython / pairwise2) |
| Medium | Sliding-window GC content plot |
| Medium | Restriction enzyme cut-site mapping |
| Medium | FASTQ support + quality score parsing |
| Low | Interactive web dashboard (Streamlit or Flask) |
| Low | RNA secondary structure prediction (ViennaRNA API) |
| Low | GenBank / EMBL flat-file parser |

---

## Dependencies

| Package | Purpose |
|---|---|
| `matplotlib` | Chart generation |
| `numpy` | Numerical support |
| `pytest` | Unit testing |

No external bioinformatics libraries are required — all core logic is implemented from scratch to demonstrate understanding of the underlying algorithms.

---

## Author

**Shubham Pradeep Pandey**  
BSc Biotechnology | Aspiring MSc Bioinformatics (Germany)  
[GitHub]: https://github.com/pp1247shu-web  [LinkedIn]: https://www.linkedin.com/in/shubham-pandey-8b544a313 

---

## License

MIT License — free to use, modify, and distribute with attribution.
