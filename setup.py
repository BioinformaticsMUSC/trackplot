"""
Setup script for the trackplot package.
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="trackplot",
    version="1.0.0",
    author="Enhanced Track Plot Generator",
    description="A comprehensive toolkit for creating publication-quality genomic track plots",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.19.0",
        "pandas>=1.3.0",
        "matplotlib>=3.3.0",
        "seaborn>=0.11.0",
        "pyBigWig>=0.3.0",
        "gffutils>=0.10.0",
    ],
    extras_require={
        "bam": ["pysam>=0.16.0"],
        "yaml": ["PyYAML>=5.4.0"],
        "dev": ["pytest>=6.0", "black", "flake8"],
    },
    entry_points={
        "console_scripts": [
            "trackplot=trackplot.cli:main",
            "trackplot-db=trackplot.db_cli:main",
        ],
    },
    include_package_data=True,
    package_data={
        'trackplot': ['data/*'],
    },
)