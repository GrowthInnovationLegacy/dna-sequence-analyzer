"""
DNA Sequence Analysis Tool - Enhanced Version
Complete implementation with all improvements integrated.
"""

import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from io import StringIO, BytesIO
import time

# Import custom modules
from src.sequence_analysis import (
    calculate_gc_content,
    calculate_sliding_window_gc,
    get_complement,
    get_reverse_complement,
    find_orfs_comprehensive,
    transcribe_dna_to_rna,
    translate_sequence,
    find_motif
)
from src.sequence_io import (
    parse_fasta,
    validate_sequence,
    SequenceRecord
)
from src.visualization import (
    plot_gc_content_window,
    plot_nucleotide_composition,
    plot_orf_locations,
    create_sequence_summary_dashboard
)
from src.utils import (
    get_sequence_statistics,
    calculate_molecular_weight,
    RESTRICTION_SITES
)

# Page configuration
st.set_page_config(
    page_title="DNA Sequence Analyzer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
    <style>
    .main-header {
        font-size: 3rem;
        font-weight: bold;
        color: #DAA520;
        text-align: center;
        padding: 1rem;
        background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
        -webkit-background-clip: text;
    }
    .metric-card {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #667eea;
    }
    .info-box {
        background-color: #FF5F15;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #3498db;
        margin: 1rem 0;
    }
    </style>
""", unsafe_allow_html=True)


def export_results_to_csv(results: dict, sequence_id: str) -> str:
    """Export analysis results to CSV format."""
    output = StringIO()
    
    # Header
    output.write(f"DNA Sequence Analysis Results\n")
    output.write(f"Sequence ID: {sequence_id}\n")
    output.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    
    # Statistics
    if 'statistics' in results:
        output.write("Sequence Statistics\n")
        for key, value in results['statistics'].items():
            output.write(f"{key},{value}\n")
        output.write("\n")
    
    # ORFs
    if 'orfs' in results and results['orfs']:
        output.write("Open Reading Frames\n")
        output.write("Frame,Start,End,Length,Protein\n")
        for orf in results['orfs']:
            output.write(f"{orf['frame']},{orf['start']},{orf['end']},{orf['length']},{orf['protein']}\n")
        output.write("\n")
    
    # Motifs
    if 'motifs' in results and results['motifs']:
        output.write("Motif Matches\n")
        output.write("Position,Sequence,Strand\n")
        for motif in results['motifs']:
            output.write(f"{motif['position']},{motif['sequence']},{motif['strand']}\n")
    
    return output.getvalue()


def main():
    """Enhanced main application function with all features."""
    
    # Header
    st.markdown('<h1 class="main-header">🧬 DNA Sequence Analyzer</h1>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="info-box">
    <b>Comprehensive bioinformatics toolkit</b> for analyzing DNA sequences with professional-grade algorithms.
    Upload FASTA files or paste sequences to get started.
    </div>
    """, unsafe_allow_html=True)
    
    # Sidebar
    st.sidebar.header("📥 Input Configuration")
    
    # Input method selection
    input_method = st.sidebar.radio(
        "Select input method:",
        ["Paste Sequence", "Upload FASTA File", "Load Example"]
    )
    
    sequences = []
    
    # Handle different input methods
    if input_method == "Paste Sequence":
        st.sidebar.markdown("**Manual Sequence Entry**")
        seq_input = st.sidebar.text_area(
            "Enter DNA sequence:",
            height=100,
            placeholder="ATGCGATCGATCG...",
            help="Paste your DNA sequence here. Only A, T, G, C characters are valid."
        )
        
        seq_id = st.sidebar.text_input(
            "Sequence ID:",
            value="my_sequence",
            help="Provide a unique identifier for this sequence"
        )
        
        if seq_input:
            # Clean and validate
            clean_seq = ''.join(seq_input.split()).upper()
            is_valid, error_msg = validate_sequence(clean_seq, 'dna')
            
            if is_valid:
                sequences = [SequenceRecord(seq_id, clean_seq)]
                st.sidebar.success(f"✅ Loaded sequence ({len(clean_seq)} bp)")
            else:
                st.sidebar.error(f"❌ {error_msg}")
    
    elif input_method == "Upload FASTA File":
        st.sidebar.markdown("**File Upload**")
        uploaded_file = st.sidebar.file_uploader(
            "Choose a FASTA file:",
            type=['fasta', 'fa', 'fna', 'txt'],
            help="Upload a file in FASTA format"
        )
        
        if uploaded_file:
            try:
                fasta_content = uploaded_file.read().decode('utf-8')
                sequences = parse_fasta(fasta_content)
                st.sidebar.success(f"✅ Loaded {len(sequences)} sequence(s)")
                
                # Show sequence selector if multiple
                if len(sequences) > 1:
                    selected_idx = st.sidebar.selectbox(
                        "Select sequence to analyze:",
                        range(len(sequences)),
                        format_func=lambda i: sequences[i].id
                    )
                    sequences = [sequences[selected_idx]]
            except Exception as e:
                st.sidebar.error(f"❌ Error: {str(e)}")
    
    else:  # Load Example
        st.sidebar.markdown("**Example Sequences**")
        example_choice = st.sidebar.selectbox(
            "Choose example:",
            ["E. coli genome fragment", "Human insulin gene", "Simple test sequence"]
        )
        
        examples = {
            "E. coli genome fragment": (
                "ecoli_fragment",
                "ATGGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGGCTAGCTAGCTAGC"
                "TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"
                "TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG"
            ),
            "Human insulin gene": (
                "insulin",
                "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAG"
                "CCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGG"
            ),
            "Simple test sequence": (
                "test",
                "ATGGCTAGCTAGCTAGCTAGCTAGCTAG"
            )
        }
        
        seq_id, seq_data = examples[example_choice]
        sequences = [SequenceRecord(seq_id, seq_data)]
        st.sidebar.success(f"✅ Loaded example: {seq_id}")
    
    # Analysis configuration
    if sequences:
        st.sidebar.markdown("---")
        st.sidebar.header("🔬 Analysis Configuration")
        
        analyses = st.sidebar.multiselect(
            "Select analyses:",
            [
                "📊 Summary Statistics",
                "🧮 GC Content Analysis",
                "🧬 Complementary Strands",
                "🔍 Open Reading Frames",
                "🎯 Motif Search",
                "📝 Transcription",
                "🧪 Translation",
                "🔬 Restriction Sites"
            ],
            default=["📊 Summary Statistics", "🧮 GC Content Analysis"]
        )
        
        # Advanced options in collapsible section
        with st.sidebar.expander("⚙️ Advanced Options"):
            show_sequence = st.checkbox("Display full sequence", value=False)
            gc_window_size = st.slider("GC window size (bp):", 10, 200, 20, 10)
            min_orf_length = st.slider("Minimum ORF length (bp):", 30, 300, 75, 15)
        
        # Main analysis section
        for seq_record in sequences:
            # Create tabs for organization
            tab1, tab2, tab3 = st.tabs(["📈 Analysis Results", "📊 Visualizations", "💾 Export Data"])
            
            with tab1:
                st.header(f"Analysis: {seq_record.id}")
                
                # Basic info
                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    st.metric("Length", f"{len(seq_record):,} bp")
                with col2:
                    gc = calculate_gc_content(seq_record.sequence)
                    st.metric("GC Content", f"{gc:.1f}%")
                with col3:
                    at = 100 - gc
                    st.metric("AT Content", f"{at:.1f}%")
                with col4:
                    mw = calculate_molecular_weight(seq_record.sequence, 'dna')
                    st.metric("MW", f"{mw/1000:.1f} kDa")
                
                # Display sequence if requested
                if show_sequence:
                    with st.expander("🔤 View Full Sequence"):
                        st.code(seq_record.sequence, language=None)
                
                # Initialize results dictionary
                results = {}
                
                # Summary Statistics
                if "📊 Summary Statistics" in analyses:
                    st.subheader("📊 Sequence Statistics")
                    stats = get_sequence_statistics(seq_record.sequence)
                    
                    col1, col2 = st.columns(2)
                    with col1:
                        st.markdown("**Nucleotide Counts:**")
                        stats_df = pd.DataFrame({
                            'Nucleotide': ['A', 'T', 'G', 'C', 'N'],
                            'Count': [stats['a_count'], stats['t_count'], 
                                     stats['g_count'], stats['c_count'], stats['n_count']],
                            'Percentage': [
                                f"{stats['a_count']/stats['length']*100:.2f}%",
                                f"{stats['t_count']/stats['length']*100:.2f}%",
                                f"{stats['g_count']/stats['length']*100:.2f}%",
                                f"{stats['c_count']/stats['length']*100:.2f}%",
                                f"{stats['n_count']/stats['length']*100:.2f}%"
                            ]
                        })
                        st.dataframe(stats_df, hide_index=True, use_container_width=True)
                    
                    with col2:
                        st.markdown("**Physical Properties:**")
                        props_df = pd.DataFrame({
                            'Property': ['Molecular Weight', 'Melting Temp', 
                                        'Purines', 'Pyrimidines'],
                            'Value': [
                                f"{stats['molecular_weight']:.2f} Da",
                                f"{stats['melting_temp']:.1f} °C",
                                f"{stats['purine_count']} bp",
                                f"{stats['pyrimidine_count']} bp"
                            ]
                        })
                        st.dataframe(props_df, hide_index=True, use_container_width=True)
                    
                    results['statistics'] = stats
                
                # GC Content Analysis
                if "🧮 GC Content Analysis" in analyses:
                    st.subheader("🧮 GC Content Distribution")
                    
                    if len(seq_record) >= gc_window_size:
                        with st.spinner("Calculating GC content..."):
                            window_results = calculate_sliding_window_gc(
                                seq_record.sequence,
                                gc_window_size
                            )
                            
                            fig = plot_gc_content_window(
                                window_results,
                                seq_record.id,
                                show_average=True
                            )
                            st.plotly_chart(fig, use_container_width=True)
                            
                            results['gc_window'] = window_results
                    else:
                        st.warning(f"⚠️ Sequence too short for window analysis (minimum: {gc_window_size} bp)")
                
                # Complementary Strands
                if "🧬 Complementary Strands" in analyses:
                    st.subheader("🧬 Complementary Strands")
                    
                    complement = get_complement(seq_record.sequence)
                    rev_complement = get_reverse_complement(seq_record.sequence)
                    
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        st.markdown("**Complement (5'→3'):**")
                        display_len = min(100, len(complement))
                        st.code(complement[:display_len] + 
                               ("..." if len(complement) > 100 else ""), 
                               language=None)
                        
                        if st.button("📋 Copy Complement", key="copy_comp"):
                            st.code(complement, language=None)
                    
                    with col2:
                        st.markdown("**Reverse Complement (5'→3'):**")
                        st.code(rev_complement[:display_len] + 
                               ("..." if len(rev_complement) > 100 else ""),
                               language=None)
                        
                        if st.button("📋 Copy Reverse Complement", key="copy_revcomp"):
                            st.code(rev_complement, language=None)
                    
                    results['complement'] = complement
                    results['reverse_complement'] = rev_complement
                
                # Open Reading Frames
                if "🔍 Open Reading Frames" in analyses:
                    st.subheader("🔍 Open Reading Frames")
                    
                    with st.spinner("Detecting ORFs..."):
                        orfs_result = find_orfs_comprehensive(seq_record.sequence, min_length=min_orf_length)
                    
                    # Extract ORF data
                    complete_orfs = orfs_result.get("complete_orfs", [])
                    partial_orfs = orfs_result.get("partial_orfs", [])
                    orf_stats = orfs_result.get("statistics", {})
                    
                    if complete_orfs:
                        st.success(f"✅ Found {len(complete_orfs)} complete ORF(s)")
                        
                        # Create DataFrame for display
                        orf_df = pd.DataFrame(complete_orfs)
                        display_orf_df = orf_df[['frame', 'start', 'end', 'length']].copy()
                        display_orf_df['protein_preview'] = orf_df['protein'].apply(
                            lambda x: x[:30] + '...' if len(x) > 30 else x
                        )
                        
                        st.dataframe(display_orf_df, use_container_width=True, hide_index=True)
                        
                        # Show detailed view
                        with st.expander("🔬 Detailed ORF Information"):
                            if complete_orfs:  # Extra safety check
                                selected_orf = st.selectbox(
                                    "Select ORF:",
                                    range(len(complete_orfs)),
                                    format_func=lambda i: f"Frame {complete_orfs[i]['frame']} at {complete_orfs[i]['start']}-{complete_orfs[i]['end']}"
                                )
                                
                                orf_detail = complete_orfs[selected_orf]
                                st.markdown(f"**Frame:** {orf_detail['frame']}")
                                st.markdown(f"**Position:** {orf_detail['start']} - {orf_detail['end']} bp")
                                st.markdown(f"**Length:** {orf_detail['length']} bp ({orf_detail['length']//3} codons)")
                                st.markdown("**Nucleotide Sequence:**")
                                st.code(orf_detail['sequence'], language=None)
                                st.markdown("**Protein Sequence:**")
                                st.code(orf_detail['protein'], language=None)
                        
                        results['orfs'] = complete_orfs
                    else:
                        st.info("ℹ️ No ORFs found matching the criteria")
                    
                    # Show partial ORFs and statistics
                    if partial_orfs:
                        st.caption(f"Detected {len(partial_orfs)} partial ORFs (sequence end fragments)")
                    if orf_stats:
                        with st.expander("📊 ORF Statistics"):
                            st.json(orf_stats)
                
                # Motif Search
                if "🎯 Motif Search" in analyses:
                    st.subheader("🎯 Motif Search")
                    
                    col1, col2 = st.columns([2, 1])
                    with col1:
                        motif_input = st.text_input(
                            "Enter motif sequence:",
                            value="ATG",
                            help="Search for specific DNA patterns"
                        )
                    with col2:
                        search_type = st.radio(
                            "Search type:",
                            ["Custom", "Restriction Sites"],
                            horizontal=True
                        )
                    
                    if search_type == "Restriction Sites":
                        enzyme = st.selectbox(
                            "Select enzyme:",
                            list(RESTRICTION_SITES.keys())
                        )
                        motif_input = RESTRICTION_SITES[enzyme]
                        st.info(f"🔬 {enzyme} recognition site: {motif_input}")
                    
                    if motif_input:
                        matches = find_motif(seq_record.sequence, motif_input)
                        
                        if matches:
                            st.success(f"✅ Found {len(matches)} occurrence(s)")
                            
                            matches_df = pd.DataFrame(matches)
                            st.dataframe(matches_df, use_container_width=True, hide_index=True)
                            
                            results['motifs'] = matches
                        else:
                            st.info("ℹ️ Motif not found in sequence")
                
                # Transcription
                if "📝 Transcription" in analyses:
                    st.subheader("📝 Transcription (DNA → RNA)")
                    
                    rna_seq = transcribe_dna_to_rna(seq_record.sequence)
                    st.markdown("**mRNA Sequence (5'→3'):**")
                    display_len = min(200, len(rna_seq))
                    st.code(rna_seq[:display_len] + 
                           ("..." if len(rna_seq) > 200 else ""),
                           language=None)
                    
                    if st.button("📋 Show Full RNA Sequence"):
                        st.code(rna_seq, language=None)
                    
                    results['rna'] = rna_seq
                
                # Translation
                if "🧪 Translation" in analyses:
                    st.subheader("🧪 Translation (DNA → Protein)")
                    
                    try:
                        # Try translation
                        if len(seq_record.sequence) % 3 != 0:
                            st.warning(f"⚠️ Sequence length ({len(seq_record)}) not multiple of 3. Trimming...")
                            trim_seq = seq_record.sequence[:len(seq_record.sequence) - (len(seq_record.sequence) % 3)]
                        else:
                            trim_seq = seq_record.sequence
                        
                        protein_seq = translate_sequence(trim_seq)
                        
                        st.markdown(f"**Protein Sequence:** ({len(protein_seq)} amino acids)")
                        st.code(protein_seq, language=None)
                        
                        # Calculate molecular weight
                        protein_mw = calculate_molecular_weight(protein_seq, 'protein')
                        st.metric("Protein Molecular Weight", f"{protein_mw/1000:.2f} kDa")
                        
                        results['protein'] = protein_seq
                        
                    except ValueError as e:
                        st.error(f"❌ Translation error: {str(e)}")
                
                # Restriction Sites
                if "🔬 Restriction Sites" in analyses:
                    st.subheader("🔬 Restriction Enzyme Sites")
                    
                    st.markdown("**Common Restriction Sites Found:**")
                    
                    found_sites = []
                    for enzyme, site in RESTRICTION_SITES.items():
                        matches = find_motif(seq_record.sequence, site)
                        if matches:
                            found_sites.append({
                                'Enzyme': enzyme,
                                'Recognition Site': site,
                                'Occurrences': len(matches),
                                'Positions': ', '.join([str(m['position']) for m in matches[:5]]) + 
                                           ('...' if len(matches) > 5 else '')
                            })
                    
                    if found_sites:
                        sites_df = pd.DataFrame(found_sites)
                        st.dataframe(sites_df, use_container_width=True, hide_index=True)
                        results['restriction_sites'] = found_sites
                    else:
                        st.info("ℹ️ No common restriction sites found")
            
            with tab2:
                st.header("📊 Visualizations")
                
                # Nucleotide composition
                st.subheader("Nucleotide Composition")
                fig_comp = plot_nucleotide_composition(seq_record.sequence)
                st.plotly_chart(fig_comp, use_container_width=True)
                
                # ORF locations if available
                if 'orfs' in results and results['orfs']:
                    st.subheader("ORF Locations")
                    fig_orf = plot_orf_locations(results['orfs'], len(seq_record))
                    st.plotly_chart(fig_orf, use_container_width=True)
                
                # Comprehensive dashboard
                st.subheader("Summary Dashboard")
                fig_dashboard = create_sequence_summary_dashboard(
                    seq_record.sequence,
                    seq_record.id
                )
                st.plotly_chart(fig_dashboard, use_container_width=True)
            
            with tab3:
                st.header("💾 Export Results")
                
                st.markdown("### Download Options")
                
                col1, col2 = st.columns(2)
                
                with col1:
                    # CSV Export
                    csv_data = export_results_to_csv(results, seq_record.id)
                    st.download_button(
                        label="📥 Download Results (CSV)",
                        data=csv_data,
                        file_name=f"{seq_record.id}_analysis.csv",
                        mime="text/csv",
                        use_container_width=True
                    )
                
                with col2:
                    # FASTA Export
                    fasta_data = seq_record.to_fasta()
                    st.download_button(
                        label="📥 Download Sequence (FASTA)",
                        data=fasta_data,
                        file_name=f"{seq_record.id}.fasta",
                        mime="text/plain",
                        use_container_width=True
                    )
                
                # Results preview
                with st.expander("📄 Preview Export Data"):
                    st.text(csv_data)
    
    else:
        # Welcome screen
        st.info("👈 Please select an input method from the sidebar to begin analysis")
        
        st.markdown("### 🚀 Quick Start Guide")
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.markdown("""
            **1. Input Data**
            - Paste sequence
            - Upload FASTA file
            - Load example
            """)
        
        with col2:
            st.markdown("""
            **2. Configure**
            - Select analyses
            - Adjust parameters
            - Set options
            """)
        
        with col3:
            st.markdown("""
            **3. Analyze**
            - View results
            - Explore visualizations
            - Export data
            """)
        
        st.markdown("### 📚 Supported Analyses")
        st.markdown("""
        - **GC Content**: Calculate nucleotide composition and distribution
        - **ORF Detection**: Find potential protein-coding regions
        - **Motif Search**: Locate specific patterns and regulatory elements
        - **Translation**: Convert DNA to protein sequences
        - **Restriction Sites**: Identify enzyme recognition sequences
        """)
    
    # Footer - always visible at bottom of sidebar
    st.sidebar.markdown("---")
    st.sidebar.markdown("""
    **DNA Sequence Analyzer v1.0**  
    Built with Python, Streamlit & Biopython  
    
    📖 [Documentation](docs/user_guide.md)  
    🐛 [Report Issue](https://github.com/GIL794/dna-sequence-analyzer/issues)  
    ⭐ [Star on GitHub](https://github.com/GIL794/dna-sequence-analyzer)
    """)


if __name__ == "__main__":
    main()
