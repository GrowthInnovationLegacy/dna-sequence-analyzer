"""
Visualization functions for DNA sequence analysis results.

This module provides plotting and charting functions for visualizing
GC content, ORF locations, motif distributions, and other analysis results.
"""

import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import pandas as pd
from typing import List, Dict
import numpy as np


def plot_gc_content_window(window_results: List[Dict], 
                           sequence_name: str = "Sequence",
                           show_average: bool = True) -> go.Figure:
    """
    Create interactive plot of GC content along sequence.
    
    Parameters
    ----------
    window_results : list of dict
        Results from calculate_sliding_window_gc.
    sequence_name : str, default='Sequence'
        Name to display in title.
    show_average : bool, default=True
        Show horizontal line at average GC content.
    
    Returns
    -------
    plotly.graph_objects.Figure
        Interactive plot object.
    
    Examples
    --------
    >>> from src.sequence_analysis import calculate_sliding_window_gc
    >>> results = calculate_sliding_window_gc("ATGC" * 100, 20)
    >>> fig = plot_gc_content_window(results)
    >>> fig.show()
    """
    positions = [r['position'] for r in window_results]
    gc_values = [r['gc_content'] for r in window_results]
    
    fig = go.Figure()
    
    # Main GC content line
    fig.add_trace(go.Scatter(
        x=positions,
        y=gc_values,
        mode='lines',
        name='GC Content',
        line=dict(color='#3498db', width=2),
        hovertemplate='Position: %{x}<br>GC Content: %{y:.2f}%<extra></extra>'
    ))
    
    # Average line
    if show_average:
        avg_gc = sum(gc_values) / len(gc_values)
        fig.add_hline(
            y=avg_gc,
            line_dash="dash",
            line_color="red",
            annotation_text=f"Average: {avg_gc:.2f}%",
            annotation_position="right"
        )
    
    fig.update_layout(
        title=f"GC Content Along {sequence_name}",
        xaxis_title="Position (bp)",
        yaxis_title="GC Content (%)",
        template="plotly_white",
        hovermode='x unified',
        height=400
    )
    
    return fig


def plot_nucleotide_composition(sequence: str) -> go.Figure:
    """
    Create bar chart of nucleotide composition.
    
    Parameters
    ----------
    sequence : str
        DNA sequence to analyze.
    
    Returns
    -------
    plotly.graph_objects.Figure
        Bar chart showing nucleotide frequencies.
    """
    sequence = sequence.upper()
    
    nucleotides = ['A', 'T', 'G', 'C']
    counts = [sequence.count(nt) for nt in nucleotides]
    percentages = [(count / len(sequence) * 100) for count in counts]
    
    colors = ['#e74c3c', '#3498db', '#f39c12', '#2ecc71']
    
    fig = go.Figure(data=[
        go.Bar(
            x=nucleotides,
            y=percentages,
            text=[f'{p:.1f}%' for p in percentages],
            textposition='outside',
            marker_color=colors,
            hovertemplate='%{x}: %{y:.2f}%<br>Count: ' + 
                         '<br>'.join([f'{c}' for c in counts]) + '<extra></extra>'
        )
    ])
    
    fig.update_layout(
        title="Nucleotide Composition",
        xaxis_title="Nucleotide",
        yaxis_title="Percentage (%)",
        template="plotly_white",
        showlegend=False,
        height=400
    )
    
    return fig


def plot_orf_locations(orfs: List[Dict], sequence_length: int) -> go.Figure:
    """
    Visualize ORF locations on sequence.
    
    Parameters
    ----------
    orfs : list of dict
        ORF results from find_orfs function.
    sequence_length : int
        Total length of the sequence.
    
    Returns
    -------
    plotly.graph_objects.Figure
        Gantt-style chart showing ORF positions.
    """
    if not orfs:
        # Return empty figure with message
        fig = go.Figure()
        fig.add_annotation(
            text="No ORFs found",
            xref="paper", yref="paper",
            x=0.5, y=0.5,
            showarrow=False,
            font=dict(size=16)
        )
        return fig
    
    # Prepare data for plotting
    frames = ['+1', '+2', '+3', '-1', '-2', '-3']
    colors = {
        '+1': '#e74c3c', '+2': '#3498db', '+3': '#2ecc71',
        '-1': '#f39c12', '-2': '#9b59b6', '-3': '#1abc9c'
    }
    
    fig = go.Figure()
    
    for frame in frames:
        frame_orfs = [orf for orf in orfs if orf['frame'] == frame]
        
        if frame_orfs:
            for orf in frame_orfs:
                fig.add_trace(go.Scatter(
                    x=[orf['start'], orf['end']],
                    y=[frame, frame],
                    mode='lines+markers',
                    name=f"Frame {frame}",
                    line=dict(color=colors[frame], width=8),
                    marker=dict(size=10, symbol='diamond'),
                    hovertemplate=f"Frame: {frame}<br>" +
                                f"Position: {orf['start']}-{orf['end']}<br>" +
                                f"Length: {orf['length']} bp<br>" +
                                f"Protein: {orf['protein'][:20]}...<extra></extra>",
                    showlegend=False if frame_orfs.index(orf) > 0 else True
                ))
    
    fig.update_layout(
        title="Open Reading Frame Locations",
        xaxis_title="Position (bp)",
        yaxis_title="Reading Frame",
        template="plotly_white",
        height=400,
        xaxis=dict(range=[0, sequence_length]),
        yaxis=dict(categoryorder='array', categoryarray=frames)
    )
    
    return fig


def plot_codon_usage(sequence: str) -> go.Figure:
    """
    Create heatmap of codon usage.
    
    Parameters
    ----------
    sequence : str
        DNA sequence (should be in-frame coding sequence).
    
    Returns
    -------
    plotly.graph_objects.Figure
        Heatmap showing codon frequencies.
    """
    from src.utils import get_codon_usage, GENETIC_CODE
    
    # Get codon usage
    try:
        codon_counts = get_codon_usage(sequence)
    except ValueError:
        # Truncate to multiple of 3
        sequence = sequence[:len(sequence) - (len(sequence) % 3)]
        codon_counts = get_codon_usage(sequence)
    
    # Convert to RNA codons for display
    rna_codons = {k.replace('T', 'U'): v for k, v in codon_counts.items()}
    
    # Group by amino acid
    aa_groups = {}
    for codon, count in rna_codons.items():
        aa = GENETIC_CODE.get(codon, 'X')
        if aa not in aa_groups:
            aa_groups[aa] = {}
        aa_groups[aa][codon] = count
    
    # Create data for heatmap
    aa_list = sorted(aa_groups.keys())
    max_codons = max(len(codons) for codons in aa_groups.values())
    
    # Prepare matrix
    codon_labels = []
    counts_matrix = []
    
    for aa in aa_list:
        codons = sorted(aa_groups[aa].keys())
        counts = [aa_groups[aa][codon] for codon in codons]
        
        codon_labels.extend([f"{aa}<br>{codon}" for codon in codons])
        counts_matrix.append(counts)
    
    fig = px.bar(
        x=codon_labels,
        y=[sum(counts) for counts in counts_matrix for _ in counts_matrix[0]],
        title="Codon Usage Frequency",
        labels={'x': 'Codon', 'y': 'Count'}
    )
    
    fig.update_layout(
        template="plotly_white",
        height=500,
        xaxis={'categoryorder': 'total descending'}
    )
    
    return fig


def plot_amino_acid_composition(protein_sequence: str) -> go.Figure:
    """
    Visualize amino acid composition of protein sequence.
    
    Parameters
    ----------
    protein_sequence : str
        Protein sequence in single-letter code.
    
    Returns
    -------
    plotly.graph_objects.Figure
        Pie chart showing amino acid distribution.
    """
    from src.utils import AMINO_ACID_PROPERTIES
    
    # Count amino acids
    aa_counts = {}
    for aa in protein_sequence:
        if aa in AMINO_ACID_PROPERTIES:
            aa_counts[aa] = aa_counts.get(aa, 0) + 1
    
    # Sort by frequency
    sorted_aas = sorted(aa_counts.items(), key=lambda x: x[1], reverse=True)
    
    labels = [f"{aa} ({AMINO_ACID_PROPERTIES[aa]['name']})" 
              for aa, _ in sorted_aas]
    values = [count for _, count in sorted_aas]
    
    fig = go.Figure(data=[go.Pie(
        labels=labels,
        values=values,
        textinfo='label+percent',
        hovertemplate='%{label}<br>Count: %{value}<br>Percentage: %{percent}<extra></extra>'
    )])
    
    fig.update_layout(
        title="Amino Acid Composition",
        template="plotly_white",
        height=500
    )
    
    return fig


def create_sequence_summary_dashboard(sequence: str, 
                                     sequence_id: str = "Sequence") -> go.Figure:
    """
    Create comprehensive dashboard with multiple visualizations.
    
    Parameters
    ----------
    sequence : str
        DNA sequence to analyze.
    sequence_id : str
        Identifier for the sequence.
    
    Returns
    -------
    plotly.graph_objects.Figure
        Multi-panel dashboard figure.
    """
    from src.utils import get_sequence_statistics
    
    stats = get_sequence_statistics(sequence)
    
    # Create subplots
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            'Nucleotide Composition',
            'Purine vs Pyrimidine',
            'Sequence Statistics',
            'GC vs AT Content'
        ),
        specs=[
            [{'type': 'bar'}, {'type': 'pie'}],
            [{'type': 'table'}, {'type': 'bar'}]
        ]
    )
    
    # 1. Nucleotide composition bar
    nucleotides = ['A', 'T', 'G', 'C']
    counts = [stats[f'{nt.lower()}_count'] for nt in nucleotides]
    
    fig.add_trace(
        go.Bar(x=nucleotides, y=counts, marker_color=['#e74c3c', '#3498db', '#f39c12', '#2ecc71']),
        row=1, col=1
    )
    
    # 2. Purine vs Pyrimidine pie
    fig.add_trace(
        go.Pie(
            labels=['Purines (A+G)', 'Pyrimidines (C+T)'],
            values=[stats['purine_count'], stats['pyrimidine_count']]
        ),
        row=1, col=2
    )
    
    # 3. Statistics table
    fig.add_trace(
        go.Table(
            header=dict(values=['Statistic', 'Value']),
            cells=dict(values=[
                ['Length (bp)', 'GC Content', 'AT Content', 'MW (kDa)', 'Tm (°C)'],
                [stats['length'], f"{stats['gc_content']}%", f"{stats['at_content']}%",
                 f"{stats['molecular_weight']/1000:.2f}", f"{stats['melting_temp']}"]
            ])
        ),
        row=2, col=1
    )
    
    # 4. GC vs AT bar
    fig.add_trace(
        go.Bar(
            x=['GC Content', 'AT Content'],
            y=[stats['gc_content'], stats['at_content']],
            marker_color=['#2ecc71', '#e74c3c']
        ),
        row=2, col=2
    )
    
    fig.update_layout(
        title_text=f"Sequence Analysis Dashboard: {sequence_id}",
        showlegend=False,
        height=800,
        template="plotly_white"
    )
    
    return fig
