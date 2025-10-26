"""
Utility functions and constants for DNA sequence analysis.

This module provides helper functions, constants (genetic code tables,
IUPAC codes), and common operations used across the package.
"""

from typing import Dict, Set
import hashlib


def process_large_sequence_chunks(sequence: str, chunk_size: int = 10000, overlap: int = 500, analysis_func = None) -> list:
    
    """
    Process massive sequences without memory overflow.
    
    Use Cases:
    - Bacterial genomes (4-10 MB)
    - Plant genomes (100+ MB) 
    - Whole chromosomes
    
    How it works:
    1. Splits sequence into overlapping chunks
    2. Processes each chunk separately
    3. Adjusts positions to global coordinates
    4. Merges results
    """
    results = []
    seq_length = len(sequence)
    
    for start in range(0, seq_length, chunk_size - overlap):
        end = min(start + chunk_size, seq_length)
        chunk = sequence[start:end]
        
        # Process chunk
        chunk_results = analysis_func(chunk)
        
        # Adjust positions to global coordinates
        for result in chunk_results:
            if 'position' in result:
                result['position'] += start
            if 'start' in result:
                result['start'] += start
                result['end'] += start
        
        results.extend(chunk_results)
        
        # Memory management
        if len(results) > 10000:
            yield results  # Return partial results
            results = []
    
    yield results
