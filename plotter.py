"""
Main plotting class for creating enhanced genomic track plots.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import FuncFormatter
import seaborn as sns
import pyBigWig
import gffutils
from typing import Optional

from .config import TrackPlotConfig, create_custom_config
from .database import load_gene_database
from .utils import (
    load_config_from_file, 
    update_config_from_dict,
    validate_file_paths,
    load_narrowpeak_file,
    smooth_signal,
    normalize_signal,
    format_genomic_position,
    mb_formatter
)


class TrackPlotter:
    """Main class for creating enhanced genomic track plots."""
    
    def __init__(self, config: Optional[TrackPlotConfig] = None, config_file: Optional[str] = None):
        """
        Initialize TrackPlotter.
        
        Parameters:
        -----------
        config : TrackPlotConfig, optional
            Configuration object. If None, creates default config.
        config_file : str, optional
            Path to config file to load. If None, auto-detects common files.
        """
        # Set up matplotlib style
        plt.style.use('default')
        sns.set_palette("husl")
        
        # Initialize config
        if config is None:
            self.config = TrackPlotConfig()
        else:
            self.config = config
            
        # Load config file if specified or auto-detect
        self._load_config_file(config_file)
        
        # Initialize gene database
        self.db = None
        self._load_gene_database()
    
    def _load_config_file(self, config_file: Optional[str] = None):
        """Load configuration from file."""
        if config_file:
            # Load specified config file
            config_dict = load_config_from_file(config_file)
            if config_dict:
                self.config = update_config_from_dict(self.config, config_dict)
        else:
            # Auto-detect common config files
            config_files_to_try = [
                "tracks_config.yaml",
                "tracks_config.yml", 
                "tracks_config.json",
                "config.yaml",
                "config.yml",
                "config.json"
            ]
            
            for config_file in config_files_to_try:
                if os.path.exists(config_file):
                    config_dict = load_config_from_file(config_file)
                    if config_dict:
                        self.config = update_config_from_dict(self.config, config_dict)
                        break
    
    def _load_gene_database(self):
        """Load gene database using centralized database manager"""
        try:
            # Check if custom database files are specified
            if hasattr(self.config, '_custom_gtf') and self.config._custom_gtf:
                print(f"Using custom GTF: {self.config._custom_gtf}")
                self.db = load_gene_database("custom", self.config._custom_gtf)
            elif hasattr(self.config, 'gtf_file') and self.config.gtf_file:
                # Legacy support - use specified GTF file
                print(f"Using legacy GTF file: {self.config.gtf_file}")
                self.db = load_gene_database("custom", self.config.gtf_file)
            else:
                # Use genome-based database
                print(f"Loading {self.config.genome} genome database...")
                self.db = load_gene_database(self.config.genome)
                
        except Exception as e:
            print(f"Error loading gene database: {e}")
            print("Continuing without gene annotations...")
            self.db = None
    
    def get_signal(self, path: str, chrom: str, start: int, end: int, 
                   smooth: bool = False, smooth_window: int = 50) -> np.ndarray:
        """Enhanced signal reading supporting both BigWig and BAM files"""
        if not os.path.exists(path):
            print(f"Warning: File {path} not found, using zeros")
            return np.zeros(end - start)
        
        try:
            # Check file extension to determine how to read it
            if path.endswith('.bw') or path.endswith('.bigwig'):
                # Read BigWig file
                bw = pyBigWig.open(path)
                values = np.array(bw.values(chrom, start, end))
                bw.close()
                
            elif path.endswith('.bam'):
                # Read BAM file using pysam (if available)
                try:
                    import pysam
                    # Create coverage from BAM file
                    bamfile = pysam.AlignmentFile(path, "rb")
                    # Get coverage for the region
                    coverage = np.zeros(end - start)
                    for pileupcolumn in bamfile.pileup(chrom, start, end):
                        if start <= pileupcolumn.pos < end:
                            coverage[pileupcolumn.pos - start] = pileupcolumn.n
                    bamfile.close()
                    values = coverage
                    print(f"Read BAM file {path} - coverage calculated")
                except ImportError:
                    print(f"Warning: pysam not available for BAM file {path}, using zeros")
                    return np.zeros(end - start)
                except Exception as e:
                    print(f"Error reading BAM file {path}: {e}")
                    return np.zeros(end - start)
            else:
                print(f"Warning: Unsupported file format for {path}, using zeros")
                return np.zeros(end - start)
            
            # Handle NaN values
            values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
            
            # Optional smoothing
            if smooth and smooth_window > 1:
                values = smooth_signal(values, smooth_window)
                
            return values
        except Exception as e:
            print(f"Error reading {path}: {e}")
            return np.zeros(end - start)
    
    def plot_gene_models(self, ax, chrom: str, start: int, end: int):
        """Enhanced gene model plotting with better styling and organization"""
        if self.db is None:
            ax.text(0.5, 0.5, "Gene database not available", 
                    transform=ax.transAxes, ha='center', va='center')
            ax.set_xlim(start, end)
            ax.set_ylim(0, 1)
            return
        
        genes_in_region = []
        transcript_count = 0
        
        # Collect all genes in region
        try:
            for gene in self.db.region(region=(chrom, start, end), featuretype="gene"):
                gene_id = gene.id
                gene_name = gene.attributes.get("gene_name", [""])[0] or gene_id
                gene_type = gene.attributes.get("gene_type", ["unknown"])[0]
                
                # Get transcripts for this gene
                transcripts = list(self.db.children(gene, featuretype="transcript"))
                
                genes_in_region.append({
                    'gene': gene,
                    'name': gene_name,
                    'type': gene_type,
                    'transcripts': transcripts
                })
                transcript_count += len(transcripts)
        
        except Exception as e:
            print(f"Error accessing gene database: {e}")
            ax.text(0.5, 0.5, f"Error: {e}", transform=ax.transAxes, ha='center')
            return
        
        if not genes_in_region:
            ax.text(0.5, 0.5, "No genes found in region", 
                    transform=ax.transAxes, ha='center', va='center')
            ax.set_xlim(start, end)
            ax.set_ylim(0, 1)
            return
        
        # Color scheme for different gene types
        gene_type_colors = {
            'protein_coding': '#2E86AB',
            'lncRNA': '#A23B72', 
            'miRNA': '#F18F01',
            'pseudogene': '#C73E1D',
            'unknown': '#666666'
        }
        
        y_pos = 0
        gene_labels_added = set()
        
        # Plot each gene and its transcripts
        for gene_info in genes_in_region:
            gene = gene_info['gene']
            gene_name = gene_info['name']
            gene_type = gene_info['type']
            transcripts = gene_info['transcripts']
            
            color = gene_type_colors.get(gene_type, gene_type_colors['unknown'])
            
            for t_i, tx in enumerate(transcripts):
                # Plot transcript backbone
                tx_start = max(start, tx.start)
                tx_end = min(end, tx.end)
                
                ax.hlines(y_pos, tx_start, tx_end, 
                         linewidth=1.5, color=color, alpha=0.7)
                
                # Plot exons as rectangles
                try:
                    exons = list(self.db.children(tx, featuretype="exon", order_by="start"))
                    for exon in exons:
                        exon_start = max(start, exon.start)
                        exon_end = min(end, exon.end)
                        if exon_end > exon_start:
                            rect = patches.Rectangle(
                                (exon_start, y_pos - 0.15),
                                exon_end - exon_start,
                                0.3,
                                facecolor=color,
                                edgecolor='white',
                                linewidth=0.5,
                                alpha=0.9
                            )
                            ax.add_patch(rect)
                except:
                    # If exon information is not available, just show the transcript
                    pass
                
                # Add gene label (only once per gene)
                if gene_name not in gene_labels_added:
                    ax.text(tx_start - (end - start) * 0.01, y_pos + 0.25, 
                           gene_name, fontsize=9, fontweight='bold',
                           verticalalignment='bottom', color=color)
                    gene_labels_added.add(gene_name)
                
                y_pos += 0.8
        
        # Formatting
        ax.set_ylim(-0.5, y_pos + 0.5)
        ax.set_xlim(start, end)
        ax.set_yticks([])
        ax.set_title("Gene Models", fontsize=12, fontweight='bold', pad=10)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.grid(True, alpha=0.3, axis='x')
    
    def create_plot(self, show_plot: bool = True, save_plot: bool = True):
        """Create the enhanced track plot"""
        print("🧬 Creating Enhanced Track Plot")
        print("=" * 50)
        
        # Validate files
        validate_file_paths(self.config)
        
        # Filter available bigwig files
        available_bigwigs = {k: v for k, v in self.config.bigwigs.items() 
                            if os.path.exists(v)}
        
        if not available_bigwigs:
            print("Warning: No BigWig files found! Please check file paths.")
            # Could still create gene model plot
        
        print(f"Found {len(available_bigwigs)} BigWig files: {list(available_bigwigs.keys())}")
        
        # Calculate number of tracks
        num_signal_tracks = len(available_bigwigs)
        num_peak_tracks = len([f for f in self.config.peak_files.values() if os.path.exists(f)]) if self.config.show_peaks else 0
        total_tracks = num_signal_tracks + (1 if num_peak_tracks > 0 else 0) + 1  # +1 for gene models
        
        # Create figure with size validation
        fig_height = (num_signal_tracks * self.config.track_height + 
                     (1 if num_peak_tracks > 0 else 0) * 1.5 + 
                     self.config.gene_track_height)
        
        # Validate figure size to prevent matplotlib errors
        max_pixels = 32767  # Matplotlib's limit per dimension
        width_pixels = self.config.figure_width * self.config.dpi
        height_pixels = fig_height * self.config.dpi
        
        if width_pixels > max_pixels or height_pixels > max_pixels:
            print(f"⚠️  Warning: Figure too large ({width_pixels:.0f}x{height_pixels:.0f} pixels)")
            print("   Reducing DPI and figure size...")
            
            # Reduce DPI first
            if self.config.dpi > 150:
                self.config.dpi = 150
            
            # Then reduce dimensions if still too large
            while (self.config.figure_width * self.config.dpi > max_pixels):
                self.config.figure_width *= 0.8
            
            while (fig_height * self.config.dpi > max_pixels):
                self.config.track_height *= 0.8
                self.config.gene_track_height *= 0.8
                fig_height = (num_signal_tracks * self.config.track_height + 
                            (1 if num_peak_tracks > 0 else 0) * 1.2 + 
                            self.config.gene_track_height)
            
            print(f"   Adjusted to: {self.config.figure_width:.1f}x{fig_height:.1f} inches, DPI={self.config.dpi}")
        
        fig, axes = plt.subplots(
            nrows=total_tracks,
            figsize=(self.config.figure_width, fig_height),
            sharex=True,
            gridspec_kw={'height_ratios': [self.config.track_height] * num_signal_tracks + 
                         ([1.2] if num_peak_tracks > 0 else []) + [self.config.gene_track_height]}
        )
        
        # Ensure axes is always a list
        if total_tracks == 1:
            axes = [axes]
        
        row = 0
        
        # Plot BigWig signal tracks
        if available_bigwigs:
            print("Plotting signal tracks...")
            x_coords = np.arange(self.config.start, self.config.end)
            
            for label, bwpath in available_bigwigs.items():
                print(f"  Processing {label}...")
                
                # Get signal
                signal = self.get_signal(
                    bwpath, self.config.chrom, self.config.start, self.config.end,
                    smooth=self.config.smooth_signal, smooth_window=self.config.smooth_window
                )
                
                # Normalize if requested
                if self.config.normalize_tracks:
                    signal = normalize_signal(signal, method='quantile')
                
                # Get color
                color = self.config.colors.get(label, f'C{row}')
                
                # Plot signal
                axes[row].fill_between(x_coords, signal, alpha=0.7, color=color, label=label)
                axes[row].plot(x_coords, signal, color=color, linewidth=1, alpha=0.9)
                
                # Formatting
                axes[row].set_ylabel(f'{label}\\nSignal', fontsize=11, fontweight='bold')
                axes[row].grid(True, alpha=0.3)
                axes[row].spines['top'].set_visible(False)
                axes[row].spines['right'].set_visible(False)
                
                # Set y-axis limits
                max_val = np.percentile(signal[signal > 0], 95) if np.any(signal > 0) else 1
                axes[row].set_ylim(0, max_val * 1.1)
                
                row += 1
        
        # Plot peaks if requested
        if self.config.show_peaks and num_peak_tracks > 0:
            print("Plotting peaks...")
            ax_peaks = axes[row]
            
            peak_y = 0
            colors_used = []
            
            for peak_name, peak_file in self.config.peak_files.items():
                if not os.path.exists(peak_file):
                    continue
                    
                # Extract track name from peak name (e.g., "3AC_peaks" -> "3AC")
                track_name = peak_name.replace('_peaks', '')
                color = self.config.colors.get(track_name, f'C{len(colors_used)}')
                colors_used.append(color)
                
                peaks_df = load_narrowpeak_file(peak_file)
                if peaks_df.empty:
                    continue
                
                # Filter peaks in region
                region_peaks = peaks_df[
                    (peaks_df['chrom'] == self.config.chrom) &
                    (peaks_df['start'] >= self.config.start) &
                    (peaks_df['end'] <= self.config.end)
                ]
                
                # Plot peaks as rectangles
                for _, peak in region_peaks.iterrows():
                    rect = patches.Rectangle(
                        (peak['start'], peak_y),
                        peak['end'] - peak['start'],
                        0.8,
                        facecolor=color,
                        alpha=self.config.peak_alpha,
                        edgecolor=color,
                        linewidth=1
                    )
                    ax_peaks.add_patch(rect)
                
                # Add label
                ax_peaks.text(self.config.start, peak_y + 0.4, track_name, 
                             fontsize=10, fontweight='bold', color=color)
                peak_y += 1
            
            ax_peaks.set_ylim(-0.2, peak_y + 0.2)
            ax_peaks.set_ylabel('Peaks', fontsize=11, fontweight='bold')
            ax_peaks.set_yticks([])
            ax_peaks.spines['top'].set_visible(False)
            ax_peaks.spines['right'].set_visible(False)
            ax_peaks.spines['left'].set_visible(False)
            ax_peaks.grid(True, alpha=0.3, axis='x')
            
            row += 1
        
        # Plot gene models
        print("Plotting gene models...")
        self.plot_gene_models(axes[row], self.config.chrom, self.config.start, self.config.end)
        
        # Format x-axis
        axes[-1].xaxis.set_major_formatter(FuncFormatter(mb_formatter))
        axes[-1].set_xlabel(f'{self.config.chrom}: {format_genomic_position(self.config.start)} - {format_genomic_position(self.config.end)}', 
                           fontsize=12, fontweight='bold')
        
        # Overall formatting
        plt.suptitle(f'Trackplot: {self.config.chrom}:{format_genomic_position(self.config.start)}-{format_genomic_position(self.config.end)}', 
                    fontsize=12, fontweight='bold', y=0.98)  # Reduced font size
        
        # Try tight_layout with error handling
        try:
            plt.tight_layout()
            plt.subplots_adjust(top=0.95)
        except Exception as e:
            print(f"Warning: Layout adjustment failed: {e}")
            # Use basic spacing instead
            plt.subplots_adjust(top=0.92, bottom=0.08, hspace=0.3)
        
        # Save with high resolution
        if save_plot:
            output_file = f"TRACK_{self.config.chrom}_{self.config.start}_{self.config.end}.pdf"
            plt.savefig(output_file, dpi=self.config.dpi, bbox_inches='tight')
            plt.savefig(output_file.replace('.pdf', '.png'), dpi=self.config.dpi, bbox_inches='tight')
            
            print(f"\nPlot saved as: {output_file}")
            print(f"Also saved as: {output_file.replace('.pdf', '.png')}")
        
        if show_plot:
            plt.show()
        
        print("\\n✅ Track plot completed successfully!")
        return fig