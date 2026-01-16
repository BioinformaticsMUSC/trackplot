"""
Command line interface for the trackplot package.
"""

import sys
from .plotter import TrackPlotter


def main():
    """Main CLI entry point"""
    # Create plotter with auto-detected config
    plotter = TrackPlotter()
    
    # Create the plot
    fig = plotter.create_plot(show_plot=True, save_plot=True)
    
    return fig


if __name__ == "__main__":
    main()