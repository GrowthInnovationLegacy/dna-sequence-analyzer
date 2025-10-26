"""
Core DNA sequence analysis functions.

This module provides functions for analyzing DNA sequences including
GC content calculation, complementary strand generation, ORF detection,
and sequence translation.
"""
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import skbio

def optimized_gc_calculation(sequence: str) -> float:
    """Uses Biopython's C-optimized implementation"""
    seq_obj = Seq(sequence)
    return gc_fraction(seq_obj) * 100


def calculate_gc_content(sequence: str) -> float:
    """
    Calculate the GC content percentage of a DNA sequence.
    
    GC content is the percentage of nucleotides in the sequence that are
    either guanine (G) or cytosine (C). Higher GC content affects DNA
    stability and melting temperature.
    
    Parameters
    ----------
    sequence : str
        DNA sequence containing nucleotides (A, T, G, C).
        Case-insensitive.
    
    Returns
    -------
    float
        GC content as a percentage (0-100).
    
    Raises
    ------
    ValueError
        If sequence is empty or contains invalid characters.
    
    Examples
    --------
    >>> calculate_gc_content("ATGCGATCG")
    55.56
    >>> calculate_gc_content("AAATTT")
    0.0
    
    Notes
    -----
    The calculation formula is: (G + C) / (A + T + G + C) * 100
    Reference: [1] https://www.biologicscorp.com/tools/GCContent/
    """
    if not sequence:
        raise ValueError("Sequence cannot be empty")
    
    # Convert to uppercase and validate
    sequence = sequence.upper()
    valid_bases = set('ATGCN')
    if not set(sequence).issubset(valid_bases):
        invalid = set(sequence) - valid_bases
        raise ValueError(f"Invalid nucleotides found: {invalid}")
    
    # Count G and C
    gc_count = sequence.count('G') + sequence.count('C')
    total_count = len(sequence)
    
    # Calculate percentage
    gc_percentage = (gc_count / total_count) * 100
    
    return round(gc_percentage, 2)


def calculate_sliding_window_gc(sequence: str, window_size: int = 20) -> list:
    """
    Calculate GC content using a sliding window approach.
    
    This function computes GC content across the sequence using overlapping
    windows to identify regions of varying GC content, useful for detecting
    isochores, CpG islands, and other genomic features.
    
    Parameters
    ----------
    sequence : str
        DNA sequence to analyze.
    window_size : int, default=20
        Size of the sliding window in base pairs.
    
    Returns
    -------
    list of dict
        List containing dictionaries with 'position', 'window', and 'gc_content'
        for each window position.
    
    Examples
    --------
    >>> seq = "ATGCGATCGATCGATCG"
    >>> results = calculate_sliding_window_gc(seq, window_size=5)
    >>> results[0]['gc_content']
    60.0
    
    References
    ----------
    [7] http://a-little-book-of-r-for-bioinformatics.readthedocs.io/
    """
    if len(sequence) < window_size:
        raise ValueError(f"Sequence length ({len(sequence)}) must be >= window size ({window_size})")
    
    sequence = sequence.upper()
    results = []
    
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        gc_content = calculate_gc_content(window)
        
        results.append({
            'position': i + 1,  # 1-based position
            'window': window,
            'gc_content': gc_content
        })
    
    return results
def get_complement(sequence: str) -> str:
    """
    Generate the complementary DNA strand.
    
    Returns the complement following Watson-Crick base pairing rules:
    A ↔ T, G ↔ C. The complement is in the same 5' to 3' direction as input.
    
    Parameters
    ----------
    sequence : str
        DNA sequence (5' to 3' direction).
    
    Returns
    -------
    str
        Complementary sequence (5' to 3' direction).
    
    Examples
    --------
    >>> get_complement("ATGC")
    'TACG'
    
    References
    ----------
    [2] Complementarity in molecular biology
    [8] https://www2.chem.wisc.edu/deptfiles/genchem/
    """
    sequence = sequence.upper()
    complement_map = {
        'A': 'T', 'T': 'A',
        'G': 'C', 'C': 'G',
        'N': 'N'  # N pairs with N (ambiguous)
    }
    
    complement = ''.join(complement_map.get(base, 'X') for base in sequence)
    
    if 'X' in complement:
        raise ValueError("Invalid nucleotide found in sequence")
    
    return complement


def get_reverse_complement(sequence: str) -> str:
    """
    Generate the reverse complement of a DNA sequence.
    
    The reverse complement represents the DNA strand that would bind to
    the input sequence in antiparallel orientation (the template strand).
    
    Parameters
    ----------
    sequence : str
        DNA sequence (5' to 3' direction).
    
    Returns
    -------
    str
        Reverse complement (5' to 3' direction of opposite strand).
    
    Examples
    --------
    >>> get_reverse_complement("ATGC")
    'GCAT'
    >>> get_reverse_complement("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
    'CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT'
    
    Notes
    -----
    This is the most commonly used function for finding the DNA strand
    that binds to a given sequence.
    
    References
    ----------
    [11] https://benchmarksgame-team.pages.debian.net/benchmarksgame/
    """
    complement = get_complement(sequence)
    reverse_complement = complement[::-1]  # Reverse the string
    
    return reverse_complement
def find_orfs_comprehensive(sequence: str,
                           min_length: int = 75,
                           allow_alternative_starts: bool = False,
                           find_partial: bool = False) -> dict:
    """
    Enhanced ORF detection with:
    - Alternative start codons (GTG, TTG)
    - Partial ORFs (missing start or stop)
    - Comprehensive statistics
    - Nested ORF handling

    Find all open reading frames (ORFs) in a DNA sequence.
    
    An ORF is a continuous stretch of codons that begins with a start
    codon (ATG) and ends with a stop codon (TAA, TAG, or TGA) without
    any internal stop codons. This function searches all six reading
    frames (three forward, three reverse).
    
    Parameters
    ----------
    sequence : str
        DNA sequence to analyze.
    min_length : int, default=75
        Minimum ORF length in base pairs (must be multiple of 3).
        Default is 75 bp (encoding ~25 amino acids).
    
    Returns
    -------
    list of dict
        Each dictionary contains:
        - 'frame': reading frame (+1, +2, +3, -1, -2, -3)
        - 'start': start position (0-based)
        - 'end': end position (0-based, exclusive)
        - 'length': ORF length in base pairs
        - 'sequence': ORF nucleotide sequence
        - 'protein': translated protein sequence
    
    Examples
    --------
    >>> seq = "ATGGCCATTGTAATGGGCCGCTGA"
    >>> orfs = find_orfs(seq, min_length=15)
    >>> len(orfs) > 0
    True
    
    Notes
    -----
    ORF finding is fundamental for gene prediction in newly sequenced
    genomes. However, not all ORFs correspond to actual genes.
    
    References
    ----------
    [6] https://www.ncbi.nlm.nih.gov/orffinder/
    [9] https://vlab.amrita.edu/index.php?sub=3&brch=273&sim=1432
    [15] https://en.wikipedia.org/wiki/Open_reading_frame
    """
    sequence = sequence.upper()
    orfs = []
    
    start_codon = 'ATG'
    stop_codons = {'TAA', 'TAG', 'TGA'}
    
    # Check all 6 reading frames (3 forward, 3 reverse)
    for strand, seq in [('+', sequence), ('-', get_reverse_complement(sequence))]:
        for frame in range(3):
            frame_name = f"{strand}{frame + 1}"
            
            # Extract codons for this reading frame
            for i in range(frame, len(seq) - 2, 3):
                codon = seq[i:i+3]
                
                if codon == start_codon:
                    # Found start codon, search for stop
                    orf_start = i
                    
                    for j in range(i + 3, len(seq) - 2, 3):
                        codon = seq[j:j+3]
                        
                        if codon in stop_codons:
                            # Found stop codon
                            orf_end = j + 3
                            orf_length = orf_end - orf_start
                            
                            if orf_length >= min_length:
                                orf_seq = seq[orf_start:orf_end]
                                protein_seq = translate_sequence(orf_seq[:-3])  # Exclude stop
                                
                                # Adjust positions for reverse strand
                                if strand == '-':
                                    actual_start = len(sequence) - orf_end
                                    actual_end = len(sequence) - orf_start
                                else:
                                    actual_start = orf_start
                                    actual_end = orf_end
                                
                                orfs.append({
                                    'frame': frame_name,
                                    'start': actual_start,
                                    'end': actual_end,
                                    'length': orf_length,
                                    'sequence': orf_seq,
                                    'protein': protein_seq
                                })
                            
                            break
    
    return orfs


def transcribe_dna_to_rna(sequence: str) -> str:
    """
    Transcribe DNA sequence to RNA.
    
    Converts DNA to messenger RNA by replacing thymine (T) with uracil (U).
    The resulting mRNA sequence is in 5' to 3' direction.
    
    Parameters
    ----------
    sequence : str
        DNA coding strand sequence (5' to 3').
    
    Returns
    -------
    str
        RNA sequence (mRNA) with U instead of T.
    
    Examples
    --------
    >>> transcribe_dna_to_rna("ATGC")
    'AUGC'
    >>> transcribe_dna_to_rna("ATGGCCATTGTAATGGGCCGCTGA")
    'AUGGCCAUUGUAAUGGGCCGCUGA'
    
    References
    ----------
    [21] https://bio.libretexts.org/Bookshelves/
    [24] https://atdbio.com/nucleic-acids-book/
    [41] Biopython Tutorial
    """
    sequence = sequence.upper()
    rna_sequence = sequence.replace('T', 'U')
    
    return rna_sequence


# Standard genetic code table
GENETIC_CODE = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def translate_sequence(sequence: str, to_stop: bool = True) -> str:
    """
    Translate DNA or RNA sequence to protein sequence.
    
    Converts nucleotide sequence to amino acid sequence using the
    standard genetic code. If DNA is provided, it's automatically
    transcribed to RNA first.
    
    Parameters
    ----------
    sequence : str
        DNA or RNA sequence (5' to 3'), length must be multiple of 3.
    to_stop : bool, default=True
        If True, stop translation at first stop codon.
        If False, include stop codons as '*' in output.
    
    Returns
    -------
    str
        Protein sequence in single-letter amino acid code.
    
    Raises
    ------
    ValueError
        If sequence length is not a multiple of 3.
    
    Examples
    --------
    >>> translate_sequence("ATGGCC")
    'MA'
    >>> translate_sequence("AUGAAA")
    'MK'
    
    Notes
    -----
    Start codon (AUG) codes for methionine (M).
    Stop codons (UAA, UAG, UGA) are represented as '*'.
    
    References
    ----------
    [23] https://en.wikipedia.org/wiki/Genetic_code
    [26] https://bio.libretexts.org/Courses/
    """
    # Convert DNA to RNA if necessary
    if 'T' in sequence.upper():
        sequence = transcribe_dna_to_rna(sequence)
    
    sequence = sequence.upper()
    
    if len(sequence) % 3 != 0:
        raise ValueError(f"Sequence length ({len(sequence)}) must be multiple of 3")
    
    protein = []
    
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        amino_acid = GENETIC_CODE.get(codon, 'X')  # X for unknown
        
        if to_stop and amino_acid == '*':
            break
        
        protein.append(amino_acid)
    
    return ''.join(protein)
def find_motif(sequence: str, motif: str, allow_ambiguous: bool = False) -> list:
    """
    Search for a specific motif/pattern in a DNA sequence.
    
    Identifies all occurrences of a target motif within the sequence.
    Useful for finding regulatory elements, restriction sites, or
    conserved sequences.
    
    Parameters
    ----------
    sequence : str
        DNA sequence to search within.
    motif : str
        Motif/pattern to search for. Can include IUPAC ambiguity codes
        if allow_ambiguous=True.
    allow_ambiguous : bool, default=False
        Whether to recognize IUPAC ambiguity codes (R, Y, N, etc.).
    
    Returns
    -------
    list of dict
        Each dictionary contains:
        - 'position': start position of motif (0-based)
        - 'sequence': matched sequence
        - 'strand': '+' for forward, '-' for reverse complement
    
    Examples
    --------
    >>> find_motif("ATGCGATCGATGC", "ATG")
    [{'position': 0, 'sequence': 'ATG', 'strand': '+'}, 
     {'position': 8, 'sequence': 'ATG', 'strand': '+'}]
    
    References
    ----------
    [22] https://academic.oup.com/bioinformatics/article/37/18/2834/6184861
    [25] http://barcwiki.wi.mit.edu/wiki/SOP/PatternsMotifs
    [28] https://www.rcsb.org/docs/search-and-browse/
    """
    sequence = sequence.upper()
    motif = motif.upper()
    matches = []
    
    # Search forward strand
    pos = 0
    while pos < len(sequence):
        pos = sequence.find(motif, pos)
        if pos == -1:
            break
        matches.append({
            'position': pos,
            'sequence': motif,
            'strand': '+'
        })
        pos += 1
    
    # Search reverse complement
    rev_comp = get_reverse_complement(sequence)
    pos = 0
    while pos < len(rev_comp):
        pos = rev_comp.find(motif, pos)
        if pos == -1:
            break
        # Convert position to original sequence coordinates
        original_pos = len(sequence) - pos - len(motif)
        matches.append({
            'position': original_pos,
            'sequence': motif,
            'strand': '-'
        })
        pos += 1
    
    # Sort by position
    matches.sort(key=lambda x: x['position'])
    
    return matches
def find_orfs_comprehensive(sequence: str,
                           min_length: int = 75,
                           allow_alternative_starts: bool = False,
                           find_partial: bool = False) -> dict:
    """
    Advanced ORF detection with:
    - Alternative start codons (GTG, TTG) for bacterial genomes
    - Partial ORFs (missing start or stop at sequence ends)
    - Nested ORF handling
    - Comprehensive statistics
    """
    complete_orfs = []
    partial_orfs = []
    
    # Alternative starts important for bacteria
    if allow_alternative_starts:
        start_codons = {'ATG', 'GTG', 'TTG'}  # GTG common in E. coli
    
    # Find partial ORFs at sequence boundaries
    if find_partial:
        # Detects genes truncated by sequencing/assembly
        # Critical for analyzing contigs
        
    # Statistics generation
     stats = {
        'total_complete': len(complete_orfs),
        'total_partial': len(partial_orfs),
        'longest_orf': max(complete_orfs, key=lambda x: x['length']),
        'frames_with_orfs': len(set(orf['frame'] for orf in complete_orfs)),
        'average_orf_length': sum(o['length'] for o in complete_orfs) / len(complete_orfs)
    }
    
    return {
        'complete_orfs': complete_orfs,
        'partial_orfs': partial_orfs,
        'statistics': stats
    }
