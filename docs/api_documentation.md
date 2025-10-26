# DNA Sequence Analyzer - API Documentation

## Module: sequence_analysis.py

### calculate_gc_content()

calculate_gc_content(sequence: str) -> float


Calculate the GC content percentage of a DNA sequence.

**Parameters:**
- `sequence` (str): DNA sequence containing A, T, G, C

**Returns:**
- `float`: GC content as percentage (0-100)

**Raises:**
- `ValueError`: If sequence is empty or contains invalid characters

**Example:**
from src.sequence_analysis import calculate_gc_content
calculate_gc_content("ATGC")
50.0


---

### get_reverse_complement()

get_reverse_complement(sequence: str) -> str


Generate the reverse complement of a DNA sequence.

**Parameters:**
- `sequence` (str): DNA sequence (5' to 3')

**Returns:**
- `str`: Reverse complement (5' to 3' of opposite strand)

**Example:**
get_reverse_complement("ATGC")
'GCAT'


---

### find_orfs()

find_orfs(sequence: str, min_length: int = 75) -> list


Find all open reading frames in a DNA sequence.

**Parameters:**
- `sequence` (str): DNA sequence
- `min_length` (int): Minimum ORF length in bp (default: 75)

**Returns:**
- `list`: List of dictionaries containing ORF information

**Dictionary Keys:**
- `frame`: Reading frame (+1, +2, +3, -1, -2, -3)
- `start`: Start position (0-based)
- `end`: End position (0-based, exclusive)
- `length`: Length in base pairs
- `sequence`: Nucleotide sequence
- `protein`: Translated protein sequence

---

## Module: sequence_io.py

### parse_fasta()

parse_fasta(fasta_string: str) -> List[SequenceRecord]


Parse FASTA formatted sequence data.

**Parameters:**
- `fasta_string` (str): FASTA formatted text

**Returns:**
- `List[SequenceRecord]`: List of parsed sequence records

**Raises:**
- `ValueError`: If FASTA format is invalid

---

### SequenceRecord Class

Container for sequence with metadata.

**Attributes:**
- `id` (str): Sequence identifier
- `sequence` (str): Nucleotide sequence
- `description` (str): Sequence description

**Methods:**
- `__len__()`: Returns sequence length
- `__str__()`: Returns FASTA formatted string
- `to_fasta(line_length=60)`: Format with line wrapping

---

## Module: utils.py

### calculate_molecular_weight()

calculate_molecular_weight(sequence: str, seq_type: str = 'dna') -> float


Calculate molecular weight of sequence.

**Parameters:**
- `sequence` (str): Sequence string
- `seq_type` (str): Type ('dna', 'rna', 'protein')

**Returns:**
- `float`: Molecular weight in Daltons

---

## Module: visualization.py

### plot_gc_content_window()

plot_gc_content_window(
window_results: List[Dict],
sequence_name: str = "Sequence",
show_average: bool = True
) -> go.Figure


Create interactive plot of GC content.

**Parameters:**
- `window_results`: Results from calculate_sliding_window_gc()
- `sequence_name`: Display name
- `show_average`: Show average line

**Returns:**
- `plotly.graph_objects.Figure`: Interactive plot

---

For complete documentation, see source code docstrings.