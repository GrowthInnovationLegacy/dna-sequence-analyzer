"""
Sequence input/output operations.

This module handles reading and writing sequence data in various formats,
particularly FASTA format which is the standard in bioinformatics.
"""
import json
from typing import List, Dict, Tuple
from io import StringIO


class SequenceRecord:
    """
    Container for a sequence with associated metadata.
    
    Attributes
    ----------
    id : str
        Sequence identifier.
    description : str
        Sequence description/annotation.
    sequence : str
        The actual nucleotide or amino acid sequence.
    """
    
    def __init__(self, id: str, sequence: str, description: str = ""):
        self.id = id
        self.sequence = sequence.upper()
        self.description = description
    
    def __len__(self):
        return len(self.sequence)
    
    def __str__(self):
        return f">{self.id} {self.description}\n{self.sequence}"
    
    def to_fasta(self, line_length: int = 60) -> str:
        """
        Format sequence as FASTA with line wrapping.
        
        Parameters
        ----------
        line_length : int, default=60
            Number of characters per line.
        
        Returns
        -------
        str
            FASTA-formatted sequence.
        """
        header = f">{self.id} {self.description}\n"
        
        # Wrap sequence
        wrapped_seq = '\n'.join(
            self.sequence[i:i+line_length] 
            for i in range(0, len(self.sequence), line_length)
        )
        
        return header + wrapped_seq


def parse_fasta(fasta_string: str) -> List[SequenceRecord]:
    """
    Parse FASTA formatted sequence data.
    
    FASTA format specifications:
    - Header line starts with '>'
    - Sequence ID is first word after '>'
    - Description is remainder of header line
    - Sequence can span multiple lines
    
    Parameters
    ----------
    fasta_string : str
        FASTA formatted text (file contents or string).
    
    Returns
    -------
    list of SequenceRecord
        Parsed sequence records.
    
    Raises
    ------
    ValueError
        If FASTA format is invalid.
    
    Examples
    --------
    >>> fasta = ">seq1 Sample sequence\\nATGC\\nGATC\\n>seq2\\nTAGC"
    >>> records = parse_fasta(fasta)
    >>> len(records)
    2
    >>> records[0].sequence
    'ATGCGATC'
    
    References
    ----------
    [41] Biopython Tutorial
    [47] https://bi1.caltech.edu/2017/code/t01_sequence_analysis.html
    """
    records = []
    current_id = None
    current_desc = ""
    current_seq = []
    
    lines = fasta_string.strip().split('\n')
    
    for line in lines:
        line = line.strip()
        
        if not line:
            continue
        
        if line.startswith('>'):
            # Save previous record if exists
            if current_id:
                records.append(
                    SequenceRecord(current_id, ''.join(current_seq), current_desc)
                )
            
            # Parse new header
            header_parts = line[1:].split(maxsplit=1)
            current_id = header_parts[0]
            current_desc = header_parts[1] if len(header_parts) > 1 else ""
            current_seq = []
        else:
            # Accumulate sequence
            current_seq.append(line)
    
    # Save last record
    if current_id:
        records.append(
            SequenceRecord(current_id, ''.join(current_seq), current_desc)
        )
    
    if not records:
        raise ValueError("No valid sequences found in FASTA format")
    
    return records

def export_results_multiformat(results: dict,
                               sequence_id: str,
                               export_format: str = 'csv') -> str:
    """
    Export to CSV, TSV, JSON, or GenBank format.
    Critical for downstream analysis integration.
    """


def validate_sequence(sequence: str, seq_type: str = 'dna') -> Tuple[bool, str]:
    """
    Validate sequence contains only valid characters.
    
    Parameters
    ----------
    sequence : str
        Sequence to validate.
    seq_type : str, default='dna'
        Type of sequence ('dna', 'rna', or 'protein').
    
    Returns
    -------
    tuple of (bool, str)
        (is_valid, error_message)
    
    Examples
    --------
    >>> validate_sequence("ATGC", "dna")
    (True, '')
    >>> validate_sequence("ATGU", "dna")
    (False, 'Invalid DNA characters: U')
    """
    sequence = sequence.upper()
    
    valid_chars = {
        'dna': set('ATGCNRYSWKMBDHV-'),
        'rna': set('AUGCNRYSWKMBDHV-'),
        'protein': set('ACDEFGHIKLMNPQRSTVWY*-')
    }
    
    allowed = valid_chars.get(seq_type.lower(), set())
    invalid = set(sequence) - allowed
    
    if invalid:
        return False, f"Invalid {seq_type.upper()} characters: {', '.join(sorted(invalid))}"
    
    return True, ""

def robust_sequence_validator(sequence: str, strict: bool = True) -> dict:
    """
    Comprehensive validation with detailed error reporting.
    Returns dict with: is_valid, cleaned_sequence, warnings, errors
    """
    warnings = []
    errors = []
    original_length = len(sequence)
    
    # Multi-layer validation
    if not strict:
        cleaned = ''.join(c for c in sequence.upper() if c in 'ATGCNRYSWKMBDHV')
        if len(cleaned) != original_length:
            removed = original_length - len(cleaned)
            warnings.append(f"Removed {removed} non-nucleotide characters")
    else:
        cleaned = sequence.upper()
    
    # Detailed error reporting
    valid_bases = set('ATGCN')
    invalid_chars = set(cleaned) - valid_bases
    
    if invalid_chars:
        errors.append(f"Invalid DNA characters: {', '.join(sorted(invalid_chars))}")
        return {'is_valid': False, 'cleaned_sequence': cleaned, 
                'warnings': warnings, 'errors': errors}
    
    # Check for ambiguous bases
    ambiguous_count = cleaned.count('N')
    if ambiguous_count > 0:
        percentage = (ambiguous_count / len(cleaned)) * 100
        warnings.append(f"Contains {ambiguous_count} ambiguous bases ({percentage:.1f}%)")
    
    return {'is_valid': True, 'cleaned_sequence': cleaned, 
            'warnings': warnings, 'errors': errors}

