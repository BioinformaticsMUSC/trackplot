"""
Enhanced ChIP-seq/CUT&Tag Track Plotting Package

A comprehensive toolkit for creating publication-quality genomic track plots
with BigWig signals, peaks, and gene annotations.
"""

from .config import TrackPlotConfig, create_custom_config
from .plotter import TrackPlotter
from .database import GeneDatabase, load_gene_database, get_db_manager
from .utils import (
    load_config_from_file,
    validate_file_paths
)

__version__ = "1.0.0"
__author__ = "Enhanced Track Plot Generator"

# Main convenience functions for easy importing
def create_trackplot(config=None, config_file=None, save_plot=True, show_plot=True):
    """
    Create a track plot with automatic configuration loading.
    
    Parameters:
    -----------
    config : TrackPlotConfig, optional
        Configuration object. If None, uses default config.
    config_file : str, optional  
        Path to config file. If None, auto-detects common config files.
    save_plot : bool, default True
        Whether to save the plot to file
    show_plot : bool, default True
        Whether to display the plot
        
    Returns:
    --------
    matplotlib.figure.Figure
        The created figure object
    """
    plotter = TrackPlotter(config=config, config_file=config_file)
    return plotter.create_plot(save_plot=save_plot, show_plot=show_plot)

def quick_plot(chrom, start, end, save_plot=True, show_plot=True, **kwargs):
    """
    Quick plot for a genomic region with minimal setup.
    
    Parameters:
    -----------
    chrom : str
        Chromosome name (e.g., 'chr1')
    start : int
        Start position
    end : int
        End position
    save_plot : bool, default True
        Whether to save the plot to file
    show_plot : bool, default True
        Whether to display the plot
    **kwargs : dict
        Additional config parameters (bigwigs, peak_files, colors, etc.)
        
    Returns:
    --------
    matplotlib.figure.Figure
        The created figure object
    """
    config = create_custom_config(chrom, start, end)
    
    # Handle reasonable defaults for figure size to prevent massive images
    if 'figure_width' not in kwargs:
        # Scale figure width based on region size, with reasonable limits
        region_size = end - start
        if region_size < 10000:  # < 10kb
            config.figure_width = 12
        elif region_size < 100000:  # < 100kb
            config.figure_width = 14
        else:  # >= 100kb
            config.figure_width = 16
    
    if 'dpi' not in kwargs:
        config.dpi = 150  # Lower default DPI to prevent huge files
    
    # Update config with any provided parameters
    for key, value in kwargs.items():
        if hasattr(config, key):
            setattr(config, key, value)
    
    # Create plotter and generate plot
    plotter = TrackPlotter(config=config)
    return plotter.create_plot(save_plot=save_plot, show_plot=show_plot)

# Database management convenience functions
def setup_genome_database(genome: str = "hg38", force: bool = False):
    """
    Download and set up genome database.
    
    Parameters:
    -----------
    genome : str, default "hg38"
        Genome assembly (hg38, hg19, mm10, mm39)
    force : bool, default False
        Force re-download and re-creation
        
    Returns:
    --------
    Path
        Path to created database
    """
    db_manager = get_db_manager()
    return db_manager.create_database(genome, force=force)

def list_available_genomes():
    """List available genome databases and their status."""
    db_manager = get_db_manager()
    return db_manager.list_available_databases()

# Make key classes available at package level
__all__ = [
    'TrackPlotConfig',
    'TrackPlotter', 
    'GeneDatabase',
    'create_trackplot',
    'quick_plot',
    'load_config_from_file',
    'create_custom_config',
    'validate_file_paths',
    'load_gene_database',
    'setup_genome_database',
    'list_available_genomes'
]