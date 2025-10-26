"""
DNA Sequence Analyzer Package

A comprehensive bioinformatics toolkit for DNA sequence analysis.

Modules
-------
sequence_analysis : Core analysis functions (GC content, ORFs, translation)
sequence_io : File I/O operations and FASTA parsing
visualization : Plotting and data visualization
utils : Helper functions and constants
"""

__version__ = "1.0.0"
__author__ = "Your Name"
__email__ = "your.email@example.com"

# Import main functions for easy access
from .sequence_analysis import (
    calculate_gc_content,
    calculate_sliding_window_gc,
    get_complement,
    get_reverse_complement,
    find_orfs,
    transcribe_dna_to_rna,
    translate_sequence,
    find_motif
)

from .sequence_io import (
    parse_fasta,
    validate_sequence,
    SequenceRecord
)

__all__ = [
    'calculate_gc_content',
    'calculate_sliding_window_gc',
    'get_complement',
    'get_reverse_complement',
    'find_orfs',
    'transcribe_dna_to_rna',
    'translate_sequence',
    'find_motif',
    'parse_fasta',
    'validate_sequence',
    'SequenceRecord'
]
