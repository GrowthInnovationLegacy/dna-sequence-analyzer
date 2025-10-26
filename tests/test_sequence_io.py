"""
Unit tests for sequence I/O operations.

Tests FASTA parsing, validation, and file handling.
"""

import pytest
from src.sequence_io import (
    parse_fasta,
    validate_sequence,
    SequenceRecord
)


class TestSequenceRecord:
    """Tests for SequenceRecord class."""
    
    def test_sequence_record_creation(self):
        """Test basic SequenceRecord creation."""
        record = SequenceRecord("seq1", "ATGC", "Test sequence")
        
        assert record.id == "seq1"
        assert record.sequence == "ATGC"
        assert record.description == "Test sequence"
        assert len(record) == 4
    
    def test_sequence_record_uppercase(self):
        """Test that sequence is converted to uppercase."""
        record = SequenceRecord("seq1", "atgc")
        assert record.sequence == "ATGC"
    
    def test_sequence_record_str(self):
        """Test string representation."""
        record = SequenceRecord("seq1", "ATGC", "Test")
        result = str(record)
        
        assert result.startswith(">seq1")
        assert "ATGC" in result
    
    def test_sequence_record_to_fasta(self):
        """Test FASTA formatting with line wrapping."""
        long_seq = "A" * 100
        record = SequenceRecord("seq1", long_seq)
        fasta = record.to_fasta(line_length=60)
        
        lines = fasta.split('\n')
        assert lines[0].startswith(">")
        assert len(lines[1]) == 60  # First sequence line
        assert len(lines[-1]) == 40  # Last sequence line (100 % 60)


class TestFASTAParsing:
    """Tests for FASTA file parsing."""
    
    def test_parse_single_sequence(self):
        """Test parsing single FASTA sequence."""
        fasta = ">seq1 Test sequence\nATGC\nGATC"
        records = parse_fasta(fasta)
        
        assert len(records) == 1
        assert records[0].id == "seq1"
        assert records[0].sequence == "ATGCGATC"
        assert records[0].description == "Test sequence"
    
    def test_parse_multiple_sequences(self):
        """Test parsing multiple FASTA sequences."""
        fasta = ">seq1\nATGC\n>seq2\nGATC"
        records = parse_fasta(fasta)
        
        assert len(records) == 2
        assert records[0].id == "seq1"
        assert records[1].id == "seq2"
    
    def test_parse_with_blank_lines(self):
        """Test parsing with blank lines."""
        fasta = ">seq1\nATGC\n\nGATC\n\n>seq2\nCCCC"
        records = parse_fasta(fasta)
        
        assert len(records) == 2
        assert records[0].sequence == "ATGCGATC"
    
    def test_parse_no_description(self):
        """Test parsing sequence without description."""
        fasta = ">seq1\nATGC"
        records = parse_fasta(fasta)
        
        assert records[0].id == "seq1"
        assert records[0].description == ""
    
    def test_parse_empty_raises_error(self):
        """Test that empty input raises error."""
        with pytest.raises(ValueError, match="No valid sequences"):
            parse_fasta("")
    
    def test_parse_no_header_raises_error(self):
        """Test that missing header raises error."""
        with pytest.raises(ValueError):
            parse_fasta("ATGC")


class TestSequenceValidation:
    """Tests for sequence validation."""
    
    def test_validate_dna_valid(self):
        """Test validation of valid DNA sequence."""
        is_valid, msg = validate_sequence("ATGC", "dna")
        assert is_valid is True
        assert msg == ""
    
    def test_validate_dna_with_n(self):
        """Test that N is allowed in DNA."""
        is_valid, msg = validate_sequence("ATGCN", "dna")
        assert is_valid is True
    
    def test_validate_dna_invalid(self):
        """Test validation catches invalid characters."""
        is_valid, msg = validate_sequence("ATGCX", "dna")
        assert is_valid is False
        assert "Invalid" in msg
        assert "X" in msg
    
    def test_validate_rna_valid(self):
        """Test validation of valid RNA sequence."""
        is_valid, msg = validate_sequence("AUGC", "rna")
        assert is_valid is True
    
    def test_validate_rna_no_thymine(self):
        """Test that T is invalid in RNA."""
        is_valid, msg = validate_sequence("ATGC", "rna")
        assert is_valid is False
        assert "T" in msg
    
    def test_validate_protein_valid(self):
        """Test validation of valid protein sequence."""
        is_valid, msg = validate_sequence("ACDEFGHIKLM", "protein")
        assert is_valid is True
    
    def test_validate_case_insensitive(self):
        """Test that validation is case-insensitive."""
        is_valid1, _ = validate_sequence("atgc", "dna")
        is_valid2, _ = validate_sequence("ATGC", "dna")
        assert is_valid1 == is_valid2 == True


class TestEdgeCases:
    """Tests for edge cases and error handling."""
    
    def test_empty_sequence_id(self):
        """Test handling of empty sequence ID."""
        record = SequenceRecord("", "ATGC")
        assert record.id == ""
    
    def test_very_long_sequence(self):
        """Test handling of long sequence."""
        long_seq = "ATGC" * 1000
        record = SequenceRecord("seq1", long_seq)
        assert len(record) == 4000
    
    def test_whitespace_in_fasta(self):
        """Test that whitespace in sequence is handled."""
        fasta = ">seq1\nATGC   \n  GATC  "
        records = parse_fasta(fasta)
        assert records[0].sequence == "ATGC   GATC  "


# Run with: pytest tests/test_sequence_io.py -v
