"""
Utility functions for file handling, validation, and data processing.
"""

import os
import json
import pandas as pd
import numpy as np
import warnings
from typing import Dict, Any, Optional
from .config import TrackPlotConfig

# Try to import yaml for config support
try:
    import yaml
    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False


def load_config_from_file(config_path: str) -> Optional[Dict[str, Any]]:
    """Load configuration from JSON or YAML file"""
    if not os.path.exists(config_path):
        print(f"Warning: Config file {config_path} not found")
        return None
    
    try:
        with open(config_path, 'r') as f:
            if config_path.endswith('.yaml') or config_path.endswith('.yml'):
                if not YAML_AVAILABLE:
                    print("Warning: PyYAML not installed, cannot read YAML config")
                    return None
                return yaml.safe_load(f)
            elif config_path.endswith('.json'):
                return json.load(f)
            else:
                print(f"Warning: Unsupported config file format: {config_path}")
                return None
    except Exception as e:
        print(f"Error loading config file {config_path}: {e}")
        return None


def update_config_from_dict(config_obj: TrackPlotConfig, config_dict: Dict[str, Any]) -> TrackPlotConfig:
    """Update config object with values from dictionary"""
    for key, value in config_dict.items():
        if hasattr(config_obj, key):
            setattr(config_obj, key, value)
        else:
            print(f"Warning: Unknown config parameter '{key}' ignored")
    
    # Update coordinates if region parameters were changed
    if any(key in config_dict for key in ['start', 'end', 'bp_shift']):
        config_obj._update_coordinates()
    
    return config_obj


def validate_file_paths(config_obj: TrackPlotConfig) -> None:
    """Validate file paths and show warnings for missing files (non-interactive)"""
    print("\n=== FILE VALIDATION ===")
    
    missing_files = []
    found_files = []
    
    # Check GTF file
    if not os.path.exists(config_obj.gtf_file):
        missing_files.append(f"GTF: {config_obj.gtf_file}")
        print(f"⚠️  GTF file not found: {config_obj.gtf_file}")
    else:
        found_files.append(f"GTF: {config_obj.gtf_file}")
    
    # Check BigWig files
    missing_bigwigs = 0
    for name, path in config_obj.bigwigs.items():
        if not os.path.exists(path):
            missing_files.append(f"BigWig {name}: {path}")
            missing_bigwigs += 1
        else:
            found_files.append(f"BigWig {name}: {path}")
    
    # Check peak files
    missing_peaks = 0
    for name, path in config_obj.peak_files.items():
        if not os.path.exists(path):
            missing_files.append(f"Peak {name}: {path}")
            missing_peaks += 1
        else:
            found_files.append(f"Peak {name}: {path}")
    
    # Summary
    if found_files:
        print(f"✅ Found {len(found_files)} files")
    
    if missing_files:
        print(f"⚠️  Missing {len(missing_files)} files:")
        for f in missing_files[:5]:  # Show first 5 missing files
            print(f"   - {f}")
        if len(missing_files) > 5:
            print(f"   ... and {len(missing_files) - 5} more")
    
    # Specific warnings
    if missing_bigwigs == len(config_obj.bigwigs):
        print("⚠️  WARNING: No BigWig files found! Tracks will be empty.")
    
    if config_obj.show_peaks and missing_peaks == len(config_obj.peak_files):
        print("⚠️  WARNING: No peak files found! Peak section will be empty.")
        print("   Consider setting config.show_peaks = False")
    
    print("=" * 25)


def load_narrowpeak_file(filepath: str) -> pd.DataFrame:
    """Load narrowPeak file and return DataFrame"""
    if not os.path.exists(filepath):
        print(f"Warning: Peak file {filepath} not found")
        return pd.DataFrame()
    
    try:
        # narrowPeak format: chrom, start, end, name, score, strand, signalValue, pValue, qValue, peak
        columns = ['chrom', 'start', 'end', 'name', 'score', 'strand', 
                  'signalValue', 'pValue', 'qValue', 'peak']
        df = pd.read_csv(filepath, sep='\t', header=None, names=columns)
        return df
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return pd.DataFrame()


def smooth_signal(signal: np.ndarray, window_size: int) -> np.ndarray:
    """Apply moving average smoothing to signal"""
    if len(signal) < window_size:
        return signal
    
    # Use pandas rolling mean for smoothing
    s = pd.Series(signal)
    smoothed = s.rolling(window=window_size, center=True, min_periods=1).mean()
    return smoothed.values


def normalize_signal(signal: np.ndarray, method: str = 'max') -> np.ndarray:
    """Normalize signal using different methods"""
    if len(signal) == 0 or np.max(signal) == 0:
        return signal
    
    if method == 'max':
        return signal / np.max(signal)
    elif method == 'zscore':
        return (signal - np.mean(signal)) / (np.std(signal) + 1e-8)
    elif method == 'quantile':
        q95 = np.percentile(signal, 95)
        return signal / (q95 + 1e-8)
    else:
        return signal


def format_genomic_position(pos: int) -> str:
    """Format genomic position with commas"""
    return f"{pos:,}"


def mb_formatter(x, pos):
    """Format axis labels in Mb"""
    return f"{x/1e6:.1f} Mb"


def configure_bam_files(config_obj: TrackPlotConfig) -> TrackPlotConfig:
    """Configure to use BAM files instead of BigWig files"""
    config_obj.bigwigs = {
        "3AC": "bams/3AC.filtered.bam",
        "OAC": "bams/OAC.filtered.bam", 
        "3K9": "bams/3K9.filtered.bam",
        "OK9": "bams/OK9.filtered.bam"
    }
    return config_obj


def configure_demo_mode(config_obj: TrackPlotConfig) -> TrackPlotConfig:
    """Configure demo mode with synthetic data"""
    config_obj.bigwigs = {"Demo_Track": "demo.bw"}  # Will create synthetic data
    config_obj.show_peaks = False  # Turn off peaks for demo
    return config_obj


def create_synthetic_data(chrom: str, start: int, end: int) -> np.ndarray:
    """Create synthetic BigWig-like data for demonstration"""
    print("Creating synthetic data for demonstration...")
    length = end - start
    
    # Create some realistic-looking ChIP-seq like data
    np.random.seed(42)  # For reproducible results
    
    # Generate base signal with some structure
    x = np.linspace(0, 4*np.pi, length)
    signal = np.abs(np.sin(x) + 0.5*np.sin(3*x) + 0.3*np.random.normal(0, 1, length))
    signal = np.maximum(signal, 0)  # No negative values
    
    # Add some peaks
    peak_positions = [length//4, length//2, 3*length//4]
    for pos in peak_positions:
        if 0 < pos < length-100:
            peak = np.exp(-((np.arange(length) - pos)**2) / (2*50**2))
            signal += peak * np.random.uniform(2, 5)
    
    return signal