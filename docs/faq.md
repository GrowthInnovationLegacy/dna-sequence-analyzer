DNA Sequence Analyzer FAQ
General
Q: What file formats are supported?
A: FASTA (.fasta, .fa, .fna, .txt). Each sequence must begin with “>id description” and be followed by one or more lines of bases.

Q: Can I analyze multiple sequences at once?
A: Yes. Upload a multi-record FASTA file. In the app, you can choose which record to analyze; the Python API supports lists of records.

Q: How large can my sequence be?
A: The app comfortably handles sequences up to a few megabases. For very large inputs, use the chunking API (process_large_sequence_chunks) to avoid memory spikes.

Functions and Outputs
Q: Why does translation stop early?
A: By default, translation stops at the first stop codon. Set to_stop=False in translate_sequence(...) to include stop codons as “*”.

Q: Why are no ORFs found?
A: Possible reasons:

No ATG start codon in-frame

ORFs are shorter than your selected minimum length

You’re on the wrong strand—try the reverse complement

For bacterial data, enable alternative starts in the comprehensive ORF finder (if you add the advanced version).

Q: What’s the difference between complement and reverse complement?
A: Complement maps bases A↔T and G↔C in the same 5'→3' direction. Reverse complement additionally reverses orientation to represent the antiparallel binding strand.

Q: Why does GC content in the plot differ slightly from the headline metric?
A: The headline metric is global GC%, while the plot uses a sliding window average. Edge windows and rounding can cause small differences.

Q: Why am I getting “length must be multiple of 3” for translation?
A: Translation requires complete codons. The app will show a warning and trim in UI; programmatically, trim or pad to a multiple of 3.

Validation and Errors
Q: “Invalid DNA characters” error appears—what should I do?
A: Clean your input (remove spaces, digits, punctuation) or use the robust validator (strict=False) to auto-clean and warn instead of failing.

Q: FASTA parser says “No valid sequences found.”
A: Ensure your file begins with “>” and that every header has at least one nucleotide line. Remove blank lines at the start of the file.

Q: Molecular weight seems off for protein sequences.
A: Ensure you pass seq_type='protein'. For ambiguous residues (“X”), the fallback uses an average residue weight. For dsDNA or circular DNA, use function parameters.

Performance and Limits
Q: The app freezes on very large sequences.
A: Use a larger GC window (e.g., 100–200), reduce motif search operations, and prefer the chunked processing API for batch analysis. Avoid loading extremely large sequences in the browser.

Q: How can I speed up GC calculations and translation?
A: Install Biopython and scikit-bio. The project already integrates Biopython for performance; keep versions in requirements.txt.

Reproducibility and Testing
Q: How do I run tests?
A: From the repository root:

bash
pytest tests/ -v
pytest tests/ --cov=src --cov-report=html
Q: How do I add my own tests?
A: Create new test files under tests/, import functions from src/, and use pytest patterns. Add sample inputs to tests/test_data/.

Data and Exports
Q: How do I export analysis results?
A: Use the Export tab in the app (CSV, FASTA). Programmatically, use export_results_multiformat(results, id, format) for CSV, JSON, GenBank, TSV.

Q: Can I search for IUPAC-ambiguous motifs?
A: The basic motif search is exact-match. Extend find_motif to interpret IUPAC codes, or convert ambiguous codes into regex before searching.

Common Pitfalls
Q: I get “NameError: get_sequence_statistics is not defined.”
A: Ensure from src.utils import get_sequence_statistics is in app.py and the function is defined in src/utils.py (see tutorial for implementation).

Q: “calculate_molecular_weight” missing or failing.
A: Restore it in src/utils.py (wrapper over Biopython with manual fallback) and import it where needed. See tutorial and prior responses for the full implementation.