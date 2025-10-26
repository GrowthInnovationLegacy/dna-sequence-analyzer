<!-- docs/tutorial.md -->
# DNA Sequence Analyzer - Complete Tutorial

## Introduction

This tutorial will guide you through using the DNA Sequence Analyzer from basic operations to advanced analysis workflows. By the end, you'll be able to analyze real genomic data confidently.

**Time Required**: 30-45 minutes  
**Prerequisites**: Basic understanding of DNA structure (A, T, G, C bases)

---

## Part 1: Getting Started (5 minutes)

### Launch the Application

Open a new terminal
Ctrl+shift+c

Create virtual environment
python -m venv venv

Activate your virtual environment
source venv/bin/activate # macOS/Linux (also GitHub CodeSpace)
venv\Scripts\activate # Windows

Start the app
streamlit run app.py

The application opens at `http://localhost:8501`

### Interface Overview

- Left Sidebar: Input configuration and analysis selection
- Main Area: Results displayed in tabs (Analysis, Visualizations, Export)
- Top Metrics: Key statistics (length, GC%, AT%, MW)

---

## Part 2: Basic Analysis (10 minutes)

### Example 1: Analyzing a Simple Sequence

Goal: Calculate GC content of a short DNA fragment

1. Input the sequence:
   - Click sidebar → "Paste Sequence"
   - Enter: `ATGGCTAGCTAGCTAGCTAG`
   - Sequence ID: `test_fragment`

2. Select analysis:
   - Check: "📊 Summary Statistics"
   - Check: "🧮 GC Content Analysis"

3. Review results:
   - Top metrics show: Length (20 bp), GC Content (~47%)
   - Statistics table shows individual base counts
   - GC content plot shows distribution

Expected Output:
- Length: 20 bp
- GC Content: ~47.5%
- A: 1, T: 10, G: 9, C: 0

---

### Example 2: Finding Complementary Strands

Goal: Generate reverse complement for primer design

1. Input: `ATGGCCATTGTAATGGGCCGCTGA`
2. Select: "🧬 Complementary Strands"
3. Result:
   - Complement (5'→3'): `TACCGGTAACATTACCCGGCGACT`
   - Reverse Complement: `TCAGCGGCCCATTACAATGGCCAT`

Use Case: The reverse complement is your primer binding sequence.

---

## Part 3: Working with FASTA Files (10 minutes)

### Example 3: Uploading Real Gene Data

Goal: Analyze human insulin gene

1. Get sample data:
   - Use included file: `data/sample_sequences/human_insulin.fasta`
   - Or download from NCBI: Search "NM_000207 insulin"

2. Upload:
   - Sidebar → "Upload FASTA File"
   - Click "Browse files"
   - Select `human_insulin.fasta`
   - See confirmation: "✅ Loaded 1 sequence(s)"

3. Run comprehensive analysis:
   - Select all available analyses
   - Adjust window size: 30 bp
   - Minimum ORF length: 100 bp

4. Interpret results:

   Summary Statistics:
   - Length: ~330 bp
   - GC Content: ~44%
   - Molecular Weight: ~102 kDa

   ORF Detection:
   - Should find 1 complete ORF in +1 frame
   - Starts at position 0 (ATG)
   - Protein sequence: "MALWMR..." (signal peptide + insulin)

   GC Content Plot:
   - Observe relatively uniform distribution
   - Slight dip in middle (AT-rich region)

---

## Part 4: Advanced Analysis (15 minutes)

### Example 4: Bacterial Genome Fragment Analysis

Goal: Analyze E. coli genome segment for genes

Using: `data/sample_sequences/ecoli_sample.fasta`

1. Upload the file

2. Configuration:
   - Minimum ORF length: 150 bp (bacterial genes are typically 300+ bp)
   - GC window: 50 bp (larger for genome analysis)

3. Key analyses:

A. ORF Detection Results:
- Expected: 2-3 ORFs
- Frame distribution: Mostly +1, +2 (coding strand)
- Longest ORF: ~300-400 bp


B. Restriction Enzyme Sites:
- Common sites found: EcoRI, BamHI
- Use for cloning strategy planning
- Note positions for primer design

C. GC Content Visualization:
- E. coli average: ~50-51% GC
- Look for peaks (coding regions) and valleys (non-coding)
- Gene-rich regions show higher GC content

4. Export results:
- Tab 3: "💾 Export Data"
- Download CSV for spreadsheet analysis
- Download JSON for Python/R scripts

---

### Example 5: Motif Searching - Finding Promoters

Goal: Locate TATA box in promoter region

Sequence: Use first 200 bp of any gene sequence

1. Select: "🎯 Motif Search"
2. Enter motif: `TATAAA` (TATA box consensus)
3. Results interpretation:
- If found at position -30 to -25 upstream of ORF start: Likely functional promoter
- If found in multiple locations: Check strand orientation
- No match: Gene may use alternative promoter (CAAT box, GC box)

Alternative motifs to search:
- `ATG`: Start codons (should match ORF starts)
- `GAATTC`: EcoRI restriction site
- `GGATCC`: BamHI restriction site
- `AATAAA`: Polyadenylation signal (eukaryotes)

---

## Part 5: Real-World Workflows (10 minutes)

### Workflow A: Gene Characterization Pipeline

Scenario: You've cloned a new gene and need to characterize it

Input sequence (from Sanger sequencing)

Run Summary Statistics → Record GC%, length

Check Complementary Strand → Verify orientation

Find ORFs → Identify coding region

Translate → Get protein sequence

Find Restriction Sites → Plan subcloning

Export all results → CSV for lab notebook



---

### Workflow B: Primer Design Support

Scenario: Design primers for PCR amplification

Input target gene sequence

Calculate GC Content → Aim for 40-60% in primer region

Get Reverse Complement → For reverse primer

Calculate Melting Temp (from statistics) → Match Tm for both primers

Check for secondary structure (manual review)

Export sequence → Send to synthesis company



Primer Design Tips:
- Primer length: 18-25 bp
- Target Tm: 55-65°C
- GC content: 40-60%
- Avoid 3' GC clamp (max 2 G/C in last 5 bases)

---

### Workflow C: Sequence Quality Control

Scenario: Verify Sanger sequencing results

Upload your sequence

Run Summary Statistics

Check for:

N bases (ambiguous calls) > 5% → Re-sequence

Unexpected length → Check for insertions/deletions

GC content far from expected → Contamination?

Compare with reference using ORF detection

Document any discrepancies



---

## Part 6: Tips & Best Practices

### Input Preparation

DO:
- Clean sequences (remove vector, adapters)
- Use correct strand orientation (5' to 3')
- Verify sequence source and date
- Save original raw data

DON'T:
- Mix DNA and RNA in analysis
- Include sequence quality scores in FASTA
- Analyze sequences with >10% N bases
- Forget to record analysis parameters

---

### Performance Optimization

For large sequences (>100 KB):
- Use larger GC window size (100+ bp)
- Increase minimum ORF length (200+ bp)
- Export results incrementally
- Consider chunking very large files

For batch analysis:
- Process one sequence at a time
- Export results after each
- Compile results in spreadsheet
- Document all parameters used

---

### Troubleshooting Common Issues

Problem: "No ORFs found"
- Solution: Lower minimum ORF length to 30 bp
- Check: Sequence orientation (try reverse complement)
- Verify: Sequence contains start codon (ATG)

Problem: GC content plot looks noisy
- Solution: Increase window size (try 50-100 bp)
- Note: Small windows are more sensitive but noisier

Problem: Translation produces stop codons (*)
- Solution: Check reading frame (+1, +2, +3)
- Verify: Sequence doesn't have frameshifts
- Note: Internal stops may indicate pseudogenes

Problem: Upload fails
- Solution: Check FASTA format (must start with >)
- Verify: File size < 200 MB
- Try: Save as plain  (.txt or .fasta)

---

## Part 7: Practice Exercises

### Exercise 1: Mystery Sequence

Analyze this sequence and answer questions:

mystery_seq
ATGGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG



Questions:
1. What is the GC content?
2. How many ORFs can you find (minimum 20 bp)?
3. What is the translated protein sequence?
4. Is this likely a real gene? Why or why not?

Answers:
1. ~47.9%
2. 1 ORF starting at position 0
3. Protein: AASSSSSSSSSSSS
4. Unlikely - highly repetitive, low complexity

---

### Exercise 2: Cloning Strategy

You want to clone this gene into a plasmid:

gene_to_clone
ATGGCCGATTAA



Tasks:
1. Find suitable restriction sites (that DON'T cut your gene)
2. Calculate the molecular weight
3. Design primers with appropriate Tm

Hints:
- Check EcoRI (GAATTC), BamHI (GGATCC)
- Use Molecular Weight from statistics
- Aim for Tm 55-60°C

---

## Part 8: Advanced Topics

### Custom Analysis Parameters

GC Window Size Selection:
- Small windows (10-20 bp): Detect local variations, CpG islands
- Medium windows (50-100 bp): Gene-level analysis
- Large windows (200+ bp): Chromosome-scale analysis

ORF Length Thresholds:
- 30 bp: Detect small peptides, regulatory ORFs
- 75 bp (default): Balance sensitivity and specificity
- 150+ bp: Focus on likely protein-coding genes

---

### Integration with Other Tools

Export to BLAST:
1. Export FASTA format
2. Go to https://blast.ncbi.nlm.nih.gov/
3. Upload your sequence
4. Run nucleotide BLAST (blastn)

Export to Multiple Sequence Alignment:
1. Export multiple sequences as FASTA
2. Use Clustal Omega or MUSCLE
3. Visualize with Jalview

Export to Python/R Analysis:
1. Export JSON format
2. Load in script:
import json
with open('results.json') as f:
data = json.load(f)
print(data['gc_content'])



---

## Conclusion

You've now learned:
- Basic sequence analysis (GC content, statistics)
- Advanced features (ORF detection, motif search)
- Real-world workflows (cloning, primer design)
- Best practices and troubleshooting

### Next Steps

1. Practice: Analyze sequences from your own research
2. Explore: Try all analysis combinations
3. Document: Keep analysis logs for reproducibility
4. Share: Export results for collaborators

### Additional Resources

- User Guide: `docs/user_guide.md`
- API Documentation: `docs/api_documentation.md`
- FAQ: `docs/faq.md`
- NCBI Tutorials: https://www.ncbi.nlm.nih.gov/guide/

Happy Analyzing! 🧬