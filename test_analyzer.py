"""
test_analyzer.py
----------------
Unit tests for the core bioinformatics logic in analyzer.py.
Run with:  python -m pytest tests/ -v
"""

import sys
from pathlib import Path

# Ensure src/ is on the path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import pytest
from analyzer import (
    validate_sequence,
    sanitize_sequence,
    count_bases,
    calculate_gc_content,
    calculate_at_content,
    reverse_complement,
    nucleotide_frequency,
    most_common_nucleotide,
    find_motif,
    analyze_sequence,
)


# --------------------------------------------------------------------------- #
#  Fixtures
# --------------------------------------------------------------------------- #

@pytest.fixture
def simple_seq():
    return "ATGCATGC"

@pytest.fixture
def gc_rich_seq():
    return "GCGCGCGC"

@pytest.fixture
def at_rich_seq():
    return "ATATATATAT"


# --------------------------------------------------------------------------- #
#  validate_sequence
# --------------------------------------------------------------------------- #

class TestValidateSequence:
    def test_valid_uppercase(self):
        valid, inv = validate_sequence("ATGC")
        assert valid is True and inv == set()

    def test_valid_lowercase(self):
        valid, inv = validate_sequence("atgc")
        assert valid is True

    def test_invalid_chars(self):
        valid, inv = validate_sequence("ATGCN")
        assert valid is False
        assert "N" in inv

    def test_numbers_invalid(self):
        valid, inv = validate_sequence("ATGC123")
        assert valid is False

    def test_empty_string(self):
        valid, inv = validate_sequence("")
        assert valid is True    # empty → no invalid chars; length handled elsewhere


# --------------------------------------------------------------------------- #
#  sanitize_sequence
# --------------------------------------------------------------------------- #

class TestSanitizeSequence:
    def test_removes_whitespace(self):
        assert sanitize_sequence("ATG\nCATG") == "ATGCATG"

    def test_converts_to_uppercase(self):
        assert sanitize_sequence("atgc") == "ATGC"

    def test_raises_on_empty(self):
        with pytest.raises(ValueError):
            sanitize_sequence("   ")


# --------------------------------------------------------------------------- #
#  count_bases
# --------------------------------------------------------------------------- #

class TestCountBases:
    def test_known_sequence(self, simple_seq):
        counts = count_bases(simple_seq)
        assert counts == {"A": 2, "T": 2, "G": 2, "C": 2}

    def test_all_same_base(self):
        counts = count_bases("AAAA")
        assert counts["A"] == 4 and counts["T"] == 0

    def test_single_base(self):
        assert count_bases("G")["G"] == 1


# --------------------------------------------------------------------------- #
#  GC / AT content
# --------------------------------------------------------------------------- #

class TestGCContent:
    def test_fifty_percent(self, simple_seq):
        assert calculate_gc_content(simple_seq) == 50.0

    def test_hundred_percent(self, gc_rich_seq):
        assert calculate_gc_content(gc_rich_seq) == 100.0

    def test_zero_percent(self, at_rich_seq):
        assert calculate_gc_content(at_rich_seq) == 0.0

    def test_at_is_complement_of_gc(self, simple_seq):
        gc = calculate_gc_content(simple_seq)
        at = calculate_at_content(simple_seq)
        assert round(gc + at, 6) == 100.0


# --------------------------------------------------------------------------- #
#  Reverse complement
# --------------------------------------------------------------------------- #

class TestReverseComplement:
    def test_classic_example(self):
        # 5'-ATCG-3' → complement TAGC → reverse 5'-CGAT-3'
        assert reverse_complement("ATCG") == "CGAT"

    def test_palindrome(self):
        # AATT → complement TTAA → reverse AATT (palindrome)
        assert reverse_complement("AATT") == "AATT"

    def test_self_complement(self):
        rc = reverse_complement("GCGC")
        assert rc == "GCGC"

    def test_length_preserved(self, simple_seq):
        assert len(reverse_complement(simple_seq)) == len(simple_seq)

    def test_double_rc_identity(self, simple_seq):
        """RC of RC should return the original sequence."""
        assert reverse_complement(reverse_complement(simple_seq)) == simple_seq


# --------------------------------------------------------------------------- #
#  Nucleotide frequency
# --------------------------------------------------------------------------- #

class TestNucleotideFrequency:
    def test_equal_distribution(self, simple_seq):
        freq = nucleotide_frequency(simple_seq)
        for base in "ATGC":
            assert freq[base] == pytest.approx(25.0)

    def test_sum_to_100(self):
        seq  = "AAATTTGGGCCC"
        freq = nucleotide_frequency(seq)
        assert sum(freq.values()) == pytest.approx(100.0)


# --------------------------------------------------------------------------- #
#  Most common nucleotide
# --------------------------------------------------------------------------- #

class TestMostCommonNucleotide:
    def test_clear_winner(self):
        base, count = most_common_nucleotide("AAATGC")
        assert base == "A" and count == 3

    def test_single_base(self):
        base, count = most_common_nucleotide("G")
        assert base == "G" and count == 1


# --------------------------------------------------------------------------- #
#  Motif finding
# --------------------------------------------------------------------------- #

class TestFindMotif:
    def test_single_occurrence(self):
        assert find_motif("ATGCATGC", "GC") == [4, 8]

    def test_no_occurrence(self):
        assert find_motif("AAAA", "GC") == []

    def test_overlapping(self):
        # AAAT: motif AA appears at positions 1, 2
        positions = find_motif("AAAAT", "AA")
        assert 1 in positions and 2 in positions

    def test_case_insensitive(self):
        assert find_motif("ATGC", "atg") == [1]


# --------------------------------------------------------------------------- #
#  Full pipeline
# --------------------------------------------------------------------------- #

class TestAnalyzeSequence:
    def test_result_keys(self, simple_seq):
        result = analyze_sequence("Test", simple_seq)
        expected_keys = {
            "name", "sequence", "length", "base_counts",
            "gc_content", "at_content", "nucleotide_frequency",
            "most_common", "reverse_complement"
        }
        assert expected_keys.issubset(result.keys())

    def test_length_matches(self, simple_seq):
        result = analyze_sequence("Test", simple_seq)
        assert result["length"] == len(simple_seq)

    def test_gc_at_sum(self, simple_seq):
        result = analyze_sequence("Test", simple_seq)
        assert result["gc_content"] + result["at_content"] == pytest.approx(100.0)
