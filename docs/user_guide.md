# DNA Sequence Analyzer - User Guide

## Table of Contents
1. [Getting Started](#getting-started)
2. [Input Methods](#input-methods)
3. [Analysis Features](#analysis-features)
4. [Interpreting Results](#interpreting-results)
5. [Exporting Data](#exporting-data)
6. [Tips and Best Practices](#tips-and-best-practices)
7. [Troubleshooting](#troubleshooting)

---

## Getting Started

### First Time Setup

1. **Install the application**
following the README.md instructions

2. **Activate your virtual environment**:
source venv/bin/activate # macOS/Linux
venv\Scripts\activate # Windows

3. **Launch the application**:
streamlit run app.py

4. The application will open in your web browser at `http://localhost:8501`

---

## Input Methods

### Method 1: Paste Sequence

Best for quick analysis of short sequences or testing.

1. Select **"Paste Sequence"** in the sidebar
2. Enter your DNA sequence in the text area
- Use only A, T, G, C characters
- Spaces and line breaks are automatically removed
- Case doesn't matter (automatically converted to uppercase)
3. Optionally provide a sequence ID for identification
4. Click outside the text area to confirm input

**Example Input:**
ATGGCTAGCTAGCTAGCTAGCTAGC
TAGCTAGCTAGCTAGCTAGCTAGC


### Method 2: Upload FASTA File

Best for analyzing multiple sequences or working with standard bioinformatics data.

1. Select **"Upload FASTA File"** in the sidebar
2. Click **"Browse files"**
3. Select your .fasta, .fa, .fna, or .txt file
4. The application automatically:
   - Parses all sequences in the file
   - Validates format and content
   - Displays number of sequences loaded

**FASTA Format Example:**
sequence1 Description of sequence 1
ATGGCTAGCTAGCTAGC
TAGCTAGCTAGCTAGC
sequence2 Description of sequence 2
GCTAGCTAGCTAGCTAG

---

## Analysis Features

### 1. GC Content Analysis

**What it does:**
- Calculates percentage of G and C nucleotides
- Generates sliding window analysis
- Visualizes GC content variation along sequence

**When to use:**
- Identifying CpG islands
- Assessing sequence composition
- Detecting genomic regions with different properties

**Parameters:**
- **Sliding Window Size**: Adjustable from 10-200 bp
  - Smaller windows: More detailed, noisier
  - Larger windows: Smoother trends, less detail

**Interpretation:**
- **High GC content (>60%)**: Often found in gene-rich regions, promoters
- **Low GC content (<40%)**: Common in introns, repetitive DNA
- **GC variation**: Can indicate gene boundaries, regulatory regions

### 2. Complementary Strand

**What it does:**
- Generates Watson-Crick complement
- Produces reverse complement

**When to use:**
- Designing primers
- Finding binding partners
- Understanding DNA structure

**Output:**
- **Complement (5'→3')**: Same direction, complementary bases
- **Reverse Complement (5'→3')**: Opposite strand as it would appear

### 3. Open Reading Frames (ORFs)

**What it does:**
- Identifies potential protein-coding regions
- Searches all six reading frames (±1, ±2, ±3)
- Translates ORFs to protein sequences

**Parameters:**
- **Minimum ORF Length**: Default 75 bp (~25 amino acids)
  - Increase to reduce false positives
  - Decrease to find smaller peptides

**When to use:**
- Gene prediction
- Finding coding potential
- Analyzing novel sequences

**Interpretation:**
- Long ORFs (>300 bp) often represent real genes
- Multiple overlapping ORFs may indicate sequencing errors
- ORFs in negative frames suggest reverse strand genes

### 4. Motif Search

**What it does:**
- Finds exact matches of specified patterns
- Searches both forward and reverse strands
- Reports all occurrences with positions

**When to use:**
- Finding restriction enzyme sites
- Locating regulatory elements (TATA box, etc.)
- Identifying conserved sequences

**Common Motifs:**
- **TATAAA**: TATA box (promoter element)
- **GAATTC**: EcoRI restriction site
- **GGATCC**: BamHI restriction site
- **ATG**: Start codon

### 5. Transcription (DNA→RNA)

**What it does:**
- Converts DNA to messenger RNA
- Replaces thymine (T) with uracil (U)

**When to use:**
- Understanding gene expression
- Preparing for translation analysis
- Educational purposes

### 6. Translation (DNA/RNA→Protein)

**What it does:**
- Translates nucleotide sequence to amino acids
- Uses standard genetic code
- Stops at first stop codon (*)

**Requirements:**
- Sequence length must be multiple of 3
- Should start with start codon (ATG) for complete protein

**When to use:**
- Predicting protein sequences
- Analyzing coding sequences
- Understanding gene function

---

## Interpreting Results

### GC Content Visualization

- **Blue line**: GC content at each position
- **Red dashed line**: Average GC content
- **Peaks**: GC-rich regions (potential genes, exons)
- **Valleys**: AT-rich regions (introns, intergenic)

### ORF Results Table

| Column | Meaning |
|--------|---------|
| Frame | Reading frame (+/-1/2/3) |
| Start | Start position (bp) |
| End | End position (bp) |
| Length | ORF length in base pairs |
| Protein | Amino acid sequence |

### Nucleotide Composition Chart

- **Red (A)**: Adenine
- **Blue (T)**: Thymine
- **Orange (G)**: Guanine
- **Green (C)**: Cytosine

Balanced composition (~25% each) indicates random sequence or diverse content.

---

## Exporting Data

### Download Options

1. **CSV Format**: Tabular data for Excel/spreadsheets
2. **FASTA Format**: Sequence data for other bioinformatics tools
3. **JSON Format**: Structured data for programming

### Screenshots

- Use browser's screenshot tool
- Or right-click charts → "Save image as"

---

## Tips and Best Practices

### Input Preparation

✅ **DO:**
- Use sequences from reliable databases (NCBI, Ensembl)
- Remove vector sequences before analysis
- Verify sequence orientation
- Use appropriate sequence length for analysis type

❌ **DON'T:**
- Mix DNA and RNA in same sequence
- Include non-nucleotide characters (numbers, spaces will be removed)
- Use extremely long sequences (>1 Mb) without chunking

### Analysis Selection

- **For quick checks**: GC Content + Complementary Strand
- **For gene prediction**: ORFs + GC Content
- **For primer design**: Complementary Strand + Motif Search
- **For functional analysis**: Translation + ORF Detection

### Performance Optimization

- **Large sequences**: Use sliding window size >50 bp
- **Multiple sequences**: Analyze separately for clarity
- **Long analyses**: Be patient, complex calculations take time

---

## Troubleshooting

### Common Errors

**"Invalid nucleotides found"**
- **Cause**: Non-ATGC characters in sequence
- **Solution**: Check for spaces, numbers, or special characters

**"Sequence length must be multiple of 3"**
- **Cause**: Translation requires complete codons
- **Solution**: Trim sequence to length divisible by 3

**"No ORFs found"**
- **Cause**: No start codon or very short sequence
- **Solution**: Lower minimum ORF length or check sequence orientation

**"FASTA format invalid"**
- **Cause**: Missing '>' header or incorrect formatting
- **Solution**: Ensure each sequence starts with '>ID description'

### Browser Issues

- **Charts not displaying**: Try Chrome or Firefox (latest versions)
- **Upload not working**: Check file size (<200 MB)
- **Interface frozen**: Refresh page (analysis will restart)

### Getting Help

1. Check this user guide
2. Review error messages carefully
3. Try with sample data provided
4. Contact support via GitHub issues

---

## Additional Resources

- **NCBI Nucleotide Database**: https://www.ncbi.nlm.nih.gov/nucleotide/
- **Ensembl Genome Browser**: https://www.ensembl.org/
- **FASTA Format Specification**: https://en.wikipedia.org/wiki/FASTA_format
- **Genetic Code Tables**: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

---

**Version**: 1.0.0  
**Last Updated**: October 2025  
**Questions?** Open an issue on GitHub