"""
Unit tests for sequence analysis functions.

These tests ensure the correctness of core DNA analysis functions
using pytest framework.
"""

import pytest
from src.sequence_analysis import (
    calculate_gc_content,
    get_complement,
    get_reverse_complement,
    find_orfs,
    transcribe_dna_to_rna,
    translate_sequence
)


class TestGCContent:
    """Tests for GC content calculation."""
    
    def test_basic_gc_content(self):
        """Test basic GC content calculation."""
        assert calculate_gc_content("ATGC") == 50.0
        assert calculate_gc_content("AAAA") == 0.0
        assert calculate_gc_content("GGGG") == 100.0
    
    def test_case_insensitive(self):
        """Test that calculation is case-insensitive."""
        assert calculate_gc_content("atgc") == calculate_gc_content("ATGC")
    
    def test_empty_sequence(self):
        """Test error handling for empty sequence."""
        with pytest.raises(ValueError, match="cannot be empty"):
            calculate_gc_content("")
    
    def test_invalid_characters(self):
        """Test error handling for invalid nucleotides."""
        with pytest.raises(ValueError, match="Invalid nucleotides"):
            calculate_gc_content("ATGCX")
    
    @pytest.mark.parametrize("sequence,expected", [
        ("ATGCGATCG", 55.56),
        ("ATATATATAT", 0.0),
        ("GCGCGCGC", 100.0),
        ("ATGN", 25.0),  # N is allowed
    ])
    def test_parametrized_gc(self, sequence, expected):
        """Test multiple sequences with expected results."""
        assert calculate_gc_content(sequence) == pytest.approx(expected, rel=1e-2)


class TestComplementary:
    """Tests for complement functions."""
    
    def test_complement(self):
        """Test basic complementation."""
        assert get_complement("ATGC") == "TACG"
        assert get_complement("AAAA") == "TTTT"
    
    def test_reverse_complement(self):
        """Test reverse complement generation."""
        assert get_reverse_complement("ATGC") == "GCAT"
        assert get_reverse_complement("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG") == \
               "CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT"
    
    def test_double_reverse_complement(self):
        """Test that double reverse complement returns original."""
        original = "ATGCGATCG"
        result = get_reverse_complement(get_reverse_complement(original))
        assert result == original


class TestORF:
    """Tests for ORF detection."""
    
    def test_simple_orf(self):
        """Test detection of simple ORF."""
        seq = "ATGGCCTAA"  # ATG-GCC-TAA (start-codon-stop)
        orfs = find_orfs(seq, min_length=9)
        assert len(orfs) >= 1
    
    def test_no_orf(self):
        """Test sequence with no valid ORF."""
        seq = "AAAAA"
        orfs = find_orfs(seq, min_length=9)
        assert len(orfs) == 0


class TestTranscriptionTranslation:
    """Tests for transcription and translation."""
    
    def test_transcription(self):
        """Test DNA to RNA transcription."""
        assert transcribe_dna_to_rna("ATGC") == "AUGC"
        assert transcribe_dna_to_rna("TTTT") == "UUUU"
    
    def test_translation_start_codon(self):
        """Test translation with start codon."""
        # ATG codes for Methionine (M)
        assert translate_sequence("ATG") == "M"
    
    def test_translation_with_stop(self):
        """Test translation stops at stop codon."""
        # ATG-AAA-TAA (M-K-*)
        assert translate_sequence("ATGAAATAA", to_stop=True) == "MK"
    
    def test_translation_invalid_length(self):
        """Test error for non-multiple of 3."""
        with pytest.raises(ValueError, match="multiple of 3"):
            translate_sequence("ATGC")  # Length 4


# Run tests with: pytest tests/test_sequence_analysis.py -v
