# Enhanced Track Plot Package

A comprehensive Python package for creating publication-quality genomic track plots from ChIP-seq, CUT&Tag, and ATAC-seq data.

## Features

- **Multi-format support**: BigWig, BAM files, and narrowPeak files
- **Smart configuration**: Auto-detects config files or uses sensible defaults
- **Peak integration**: Visualizes peak regions as dedicated track sections
- **Gene annotation**: Enhanced gene model display with exon/intron structure
- **Publication ready**: High-resolution output in PDF and PNG formats
- **Flexible styling**: Customizable colors, normalization, and smoothing

## Quick Start

### Installation

```bash
# Install the package in development mode
pip install -e .

# With optional dependencies
pip install -e .[bam,yaml]
```

### Basic Usage

```python
import trackplot

# Simple plot with auto-detected config
fig = trackplot.create_trackplot()

# Quick plot for a specific region
fig = trackplot.quick_plot("chr1", 1000000, 1200000, 
                          bigwigs={"Sample1": "path/to/sample1.bw"})

# Advanced usage with custom config
config = trackplot.TrackPlotConfig()
config.set_region("chr8", 22070000, 22080000)
config.colors = {"3AC": "#FF6B6B", "OAC": "#4ECDC4"}
fig = trackplot.create_trackplot(config)
```

### Command Line

```bash
# Basic usage (auto-detects config files)
python -m trackplot

# Or if installed as package
trackplot
```

## Configuration

The package automatically detects these config files:
- `tracks_config.yaml/yml/json`
- `config.yaml/yml/json`

Example config file:
```yaml
chrom: "chr8"
start: 22073163
end: 22073831
bp_shift: 10000

bigwigs:
  "3AC": "path/to/3AC.bw"
  "OAC": "path/to/OAC.bw"

peak_files:
  "3AC_peaks": "beds/3AC_peaks.narrowPeak"
  "OAC_peaks": "beds/OAC_peaks.narrowPeak"

colors:
  "3AC": "#1f77b4"
  "OAC": "#ff7f0e"

normalize_tracks: true
show_peaks: true
```

## Dependencies

### Required
- numpy >= 1.19.0
- pandas >= 1.3.0  
- matplotlib >= 3.3.0
- seaborn >= 0.11.0
- pyBigWig >= 0.3.0
- gffutils >= 0.10.0

### Optional
- pysam >= 0.16.0 (for BAM file support)
- PyYAML >= 5.4.0 (for YAML config files)