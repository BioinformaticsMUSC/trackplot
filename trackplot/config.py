"""
Configuration classes and utilities for track plotting.
"""

import os
from typing import Dict, Any, Optional


class TrackPlotConfig:
    """Configuration class for track plotting with sensible defaults."""
    
    def __init__(self, genome: str = "hg38"):
        # Genomic region
        self.chrom = "chr8"
        self.start = 22073163
        self.end = 22073831
        self.bp_shift = 10000  # bases to extend on each side
        
        # Genome and database settings
        self.genome = genome  # hg38, hg19, mm10, mm39, or "custom"
        self.gtf_file = None  # Will be set by database manager
        self.db_file = None   # Will be set by database manager
        
        # Legacy support - if these are set manually, they override genome setting
        self._custom_gtf = None
        self._custom_db = None
        
        # BigWig files - UPDATE THESE PATHS TO YOUR ACTUAL BIGWIG LOCATIONS
        # Currently looking for .bw files, but you may have .bam files or need to create BigWigs
        self.bigwigs = {
            # Option 1: If you have BigWig files somewhere else, update these paths
            "3AC": "/path/to/your/3AC.bw",
            "OAC": "/path/to/your/OAC.bw", 
            "3K9": "/path/to/your/3K9.bw",
            "OK9": "/path/to/your/OK9.bw",
            
            # Option 2: If you want to use BAM files (we'll need to convert or use different method)
            # "3AC": "bams/3AC.filtered.bam",
            # "OAC": "bams/OAC.filtered.bam",
            # "3K9": "bams/3K9.filtered.bam", 
            # "OK9": "bams/OK9.filtered.bam"
        }
        
        # Peak files (narrowPeak format)
        self.peak_files = {
            "3AC_peaks": "beds/3AC_peaks.narrowPeak",
            "OAC_peaks": "beds/OAC_peaks.narrowPeak",
            "3K9_peaks": "beds/3K9_peaks.narrowPeak",
            "OK9_peaks": "beds/OK9_peaks.narrowPeak"
        }
        
        # Visualization settings
        self.colors = {
            "3AC": "#1f77b4",  # Blue
            "OAC": "#ff7f0e",  # Orange
            "3K9": "#2ca02c",  # Green
            "OK9": "#d62728"   # Red
        }
        
        self.figure_width = 12  # Reduced from 16 to prevent huge files
        self.track_height = 2.0  # Reduced from 2.5
        self.gene_track_height = 2.5  # Reduced from 3
        self.dpi = 150  # Reduced from 300 to prevent massive images
        
        # Analysis settings
        self.normalize_tracks = True
        self.show_peaks = True
        self.peak_alpha = 0.3
        self.smooth_signal = True
        self.smooth_window = 50
        
        # Update genomic coordinates
        self._update_coordinates()
    
    def _update_coordinates(self):
        """Update genomic coordinates based on bp_shift."""
        self.start = max(0, self.start - self.bp_shift)
        self.end = self.end + self.bp_shift
    
    def set_region(self, chrom: str, start: int, end: int, bp_shift: Optional[int] = None):
        """
        Set genomic region for plotting.
        
        Parameters:
        -----------
        chrom : str
            Chromosome name
        start : int
            Start position
        end : int
            End position  
        bp_shift : int, optional
            Bases to extend on each side. If None, uses current bp_shift.
        """
        self.chrom = chrom
        original_start = start
        original_end = end
        
        if bp_shift is not None:
            self.bp_shift = bp_shift
            
        # Store original coordinates before extension
        self._original_start = original_start
        self._original_end = original_end
        
        # Apply extension
        self.start = max(0, original_start - self.bp_shift)
        self.end = original_end + self.bp_shift
    
    def set_custom_database(self, gtf_file: str, db_file: str = None):
        """
        Set custom GTF and database files.
        
        Parameters:
        -----------
        gtf_file : str
            Path to GTF file
        db_file : str, optional
            Path to database file. If None, will be auto-generated.
        """
        self._custom_gtf = gtf_file
        self._custom_db = db_file
        self.genome = "custom"
    
    def copy(self):
        """Create a copy of the configuration."""
        new_config = TrackPlotConfig(genome=self.genome)
        for attr in dir(self):
            if not attr.startswith('_') and not callable(getattr(self, attr)):
                setattr(new_config, attr, getattr(self, attr))
        return new_config


def create_custom_config(chrom: str, start: int, end: int, bp_shift: int = 5000) -> TrackPlotConfig:
    """
    Create a custom configuration for a specific region.
    
    Parameters:
    -----------
    chrom : str
        Chromosome name
    start : int
        Start position
    end : int
        End position
    bp_shift : int, default 5000
        Bases to extend on each side
        
    Returns:
    --------
    TrackPlotConfig
        Configured object for the specified region
    """
    config = TrackPlotConfig()
    config.set_region(chrom, start, end, bp_shift)
    return config