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
from Bio.SeqUtils import molecular_weight as bp_molecular_weight
from Bio.Seq import Seq

def calculate_molecular_weight(sequence: str,
                               seq_type: str = "dna",
                               double_stranded: bool = False,
                               circular: bool = False,
                               monoisotopic: bool = False) -> float:
    """
    Calculate the molecular weight (Daltons) of a DNA, RNA, or protein sequence.

    Automatically uses Biopython's optimized function if available, and falls
    back to internal constants if the sequence includes ambiguous bases.

    Parameters
    ----------
    sequence : str
        Sequence string (A, T, G, C for DNA; A, U, G, C for RNA; amino acids for protein)
    seq_type : str
        'dna', 'rna', or 'protein' (default = 'dna')
    double_stranded : bool
        If True, computes MW of double-stranded sequence (2x)
    circular : bool
        If True, excludes terminal 5' and 3' groups (for plasmids)
    monoisotopic : bool
        If True, uses monoisotopic mass instead of average mass tables

    Returns
    -------
    float
        Molecular weight in Daltons (g/mol)
    """
    # Clean input
    sequence = sequence.strip().upper()
    seq_type = seq_type.lower()

    try:
        # Use Biopython implementation if possible
        mw = bp_molecular_weight(
            Seq(sequence),
            seq_type=seq_type.upper(),
            double_stranded=double_stranded,
            circular=circular,
            monoisotopic=monoisotopic
        )
    except Exception:
        # Manual fallback for ambiguous sequences
        weights = {
            "dna": {"A": 331.0, "T": 322.0, "G": 347.0, "C": 307.0, "N": 325.0},
            "rna": {"A": 329.0, "U": 306.0, "G": 345.0, "C": 305.0, "N": 321.0},
            "protein": {
                "A": 89.09, "R": 174.20, "N": 132.12, "D": 133.10,
                "C": 121.16, "E": 147.13, "Q": 146.15, "G": 75.07,
                "H": 155.16, "I": 131.18, "L": 131.18, "K": 146.19,
                "M": 149.21, "F": 165.19, "P": 115.13, "S": 105.09,
                "T": 119.12, "W": 204.23, "Y": 181.19, "V": 117.15,
                "X": 135.00  # average for unknown residues
            }
        }

        base_weights = weights.get(seq_type, weights["dna"])
        mw = sum(base_weights.get(base, 0.0) for base in sequence)

        # Adjust for double-stranded or circular DNA
        if double_stranded:
            mw *= 2
        if circular:
            # Small correction for missing terminal groups
            mw -= 36.04  # g/mol adjustment for 5'/3' phosphate loss

    return round(mw, 2)
def get_sequence_statistics(sequence: str) -> dict:
    """
    Calculate base composition and key statistics for a DNA sequence.

    This function provides detailed information including base counts,
    GC/AT content, purine/pyrimidine counts, molecular weight, and
    estimated melting temperature.

    Parameters
    ----------
    sequence : str
        DNA sequence (expected characters: A, T, G, C, N).

    Returns
    -------
    dict
        Dictionary with calculated sequence metrics.

        Keys include:
        - length : Total sequence length (bp)
        - a_count, t_count, g_count, c_count, n_count : Nucleotide counts
        - gc_content : GC content (%)
        - at_content : AT content (%)
        - molecular_weight : Estimated molecular weight in Daltons
        - melting_temp : Estimated melting temperature (°C)
        - purine_count : Total purines (A + G)
        - pyrimidine_count : Total pyrimidines (C + T)
    """
    from src.utils import calculate_molecular_weight, calculate_melting_temperature

    sequence = sequence.strip().upper()
    if not sequence:
        raise ValueError("Empty sequence: no data to analyze")

    valid_bases = set('ATGCN')
    if not set(sequence).issubset(valid_bases):
        invalid = set(sequence) - valid_bases
        raise ValueError(f"Invalid nucleotide(s) detected: {', '.join(invalid)}")

    length = len(sequence)

    # Base counts
    a_count = sequence.count('A')
    t_count = sequence.count('T')
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    n_count = sequence.count('N')

    gc_count = g_count + c_count
    at_count = a_count + t_count

    # Percent composition
    gc_content = (gc_count / length * 100) if length > 0 else 0
    at_content = (at_count / length * 100) if length > 0 else 0

    # Derived properties
    mol_weight = calculate_molecular_weight(sequence, seq_type='dna')
    melt_temp = calculate_melting_temperature(sequence)

    purine_count = a_count + g_count
    pyrimidine_count = c_count + t_count

    return {
        "length": length,
        "a_count": a_count,
        "t_count": t_count,
        "g_count": g_count,
        "c_count": c_count,
        "n_count": n_count,
        "gc_content": round(gc_content, 2),
        "at_content": round(at_content, 2),
        "molecular_weight": round(mol_weight, 2),
        "melting_temp": round(melt_temp, 2),
        "purine_count": purine_count,
        "pyrimidine_count": pyrimidine_count
    }
def calculate_melting_temperature(sequence: str, method: str = "wallace") -> float:
    """
    Estimate the DNA melting temperature (Tm).

    Parameters
    ----------
    sequence : str
        DNA sequence (must contain A, T, G, and C).
    method : str
        Calculation method ('wallace' or 'gc_content').

    Returns
    -------
    float
        Estimated melting temperature (°C).
    """
    sequence = sequence.strip().upper()
    if not sequence:
        return 0.0
    a, t, g, c = sequence.count('A'), sequence.count('T'), sequence.count('G'), sequence.count('C')

    if method == 'wallace':
        return round(2 * (a + t) + 4 * (g + c), 2)
    elif method == 'gc_content':
        gc = g + c
        length = len(sequence)
        if length == 0:
            return 0.0
        return round(64.9 + 41 * (gc - 16.4) / length, 2)
    else:
        raise ValueError("Unknown method. Use 'wallace' or 'gc_content'.")
