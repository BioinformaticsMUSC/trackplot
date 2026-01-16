"""
Database management for gene annotations.
Handles downloading, caching, and loading gene databases.
"""

import os
import gffutils
from pathlib import Path
from typing import Optional
import urllib.request
import gzip
import shutil


class GeneDatabase:
    """Manages gene annotation databases with automatic downloading and caching."""
    
    def __init__(self, package_data_dir: Optional[str] = None):
        """
        Initialize database manager.
        
        Parameters:
        -----------
        package_data_dir : str, optional
            Directory for storing databases. If None, uses package data directory.
        """
        if package_data_dir is None:
            # Use package data directory
            self.data_dir = Path(__file__).parent / "data"
        else:
            self.data_dir = Path(package_data_dir)
        
        self.data_dir.mkdir(exist_ok=True)
        
        # Common genome database configurations
        self.genome_configs = {
            "hg38": {
                "gtf_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz",
                "gtf_file": self.data_dir / "gencode.v44.annotation.gtf",
                "db_file": self.data_dir / "gencode_hg38.db"
            },
            "hg19": {
                "gtf_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz",
                "gtf_file": self.data_dir / "gencode.v19.annotation.gtf",
                "db_file": self.data_dir / "gencode_hg19.db"
            },
            "mm10": {
                "gtf_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz",
                "gtf_file": self.data_dir / "gencode.vM25.annotation.gtf",
                "db_file": self.data_dir / "gencode_mm10.db"
            },
            "mm39": {
                "gtf_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.annotation.gtf.gz",
                "gtf_file": self.data_dir / "gencode.vM32.annotation.gtf", 
                "db_file": self.data_dir / "gencode_mm39.db"
            }
        }
    
    def download_gtf(self, genome: str, force: bool = False) -> Path:
        """
        Download GTF file for specified genome.
        
        Parameters:
        -----------
        genome : str
            Genome assembly (hg38, hg19, mm10, mm39)
        force : bool, default False
            Force re-download even if file exists
            
        Returns:
        --------
        Path
            Path to downloaded GTF file
        """
        if genome not in self.genome_configs:
            raise ValueError(f"Unsupported genome: {genome}. Supported: {list(self.genome_configs.keys())}")
        
        config = self.genome_configs[genome]
        gtf_file = config["gtf_file"]
        
        if gtf_file.exists() and not force:
            print(f"GTF file already exists: {gtf_file}")
            return gtf_file
        
        print(f"Downloading GTF for {genome}...")
        gtf_url = config["gtf_url"]
        
        try:
            # Download compressed GTF
            compressed_file = gtf_file.with_suffix('.gtf.gz')
            urllib.request.urlretrieve(gtf_url, compressed_file)
            
            # Decompress
            print("Decompressing GTF file...")
            with gzip.open(compressed_file, 'rb') as f_in:
                with open(gtf_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            
            # Remove compressed file
            compressed_file.unlink()
            
            print(f"GTF downloaded: {gtf_file}")
            return gtf_file
            
        except Exception as e:
            print(f"Error downloading GTF: {e}")
            raise
    
    def create_database(self, genome: str, gtf_file: Optional[str] = None, force: bool = False) -> Path:
        """
        Create gene database from GTF file.
        
        Parameters:
        -----------
        genome : str
            Genome assembly (hg38, hg19, mm10, mm39) or "custom"
        gtf_file : str, optional
            Custom GTF file path. Required if genome="custom"
        force : bool, default False
            Force recreation even if database exists
            
        Returns:
        --------
        Path
            Path to created database file
        """
        if genome == "custom":
            if gtf_file is None:
                raise ValueError("gtf_file must be provided for custom genome")
            db_file = self.data_dir / f"custom_{Path(gtf_file).stem}.db"
            gtf_path = Path(gtf_file)
        else:
            if genome not in self.genome_configs:
                raise ValueError(f"Unsupported genome: {genome}")
            
            config = self.genome_configs[genome]
            db_file = config["db_file"]
            gtf_path = config["gtf_file"]
            
            # Download GTF if it doesn't exist
            if not gtf_path.exists():
                gtf_path = self.download_gtf(genome)
        
        if db_file.exists() and not force:
            print(f"Database already exists: {db_file}")
            return db_file
        
        print(f"Creating gene database from {gtf_path}...")
        
        try:
            gffutils.create_db(
                str(gtf_path),
                str(db_file),
                force=True,
                merge_strategy="merge",
                disable_infer_transcripts=True,
                disable_infer_genes=True
            )
            print(f"Database created: {db_file}")
            return db_file
            
        except Exception as e:
            print(f"Error creating database: {e}")
            raise
    
    def load_database(self, genome: str, gtf_file: Optional[str] = None, auto_create: bool = True) -> Optional[gffutils.FeatureDB]:
        """
        Load gene database, creating if necessary.
        
        Parameters:
        -----------
        genome : str
            Genome assembly (hg38, hg19, mm10, mm39) or "custom"
        gtf_file : str, optional
            Custom GTF file path. Required if genome="custom"
        auto_create : bool, default True
            Automatically create database if it doesn't exist
            
        Returns:
        --------
        gffutils.FeatureDB or None
            Loaded database object, or None if loading failed
        """
        try:
            if genome == "custom":
                if gtf_file is None:
                    raise ValueError("gtf_file must be provided for custom genome")
                db_file = self.data_dir / f"custom_{Path(gtf_file).stem}.db"
            else:
                if genome not in self.genome_configs:
                    raise ValueError(f"Unsupported genome: {genome}")
                db_file = self.genome_configs[genome]["db_file"]
            
            # Create database if it doesn't exist and auto_create is True
            if not db_file.exists() and auto_create:
                db_file = self.create_database(genome, gtf_file)
            
            if not db_file.exists():
                print(f"Database file not found: {db_file}")
                return None
            
            print(f"Loading gene database: {db_file}")
            db = gffutils.FeatureDB(str(db_file))
            return db
            
        except Exception as e:
            print(f"Error loading database: {e}")
            return None
    
    def list_available_databases(self) -> dict:
        """
        List available genome databases and their status.
        
        Returns:
        --------
        dict
            Dictionary with genome names and their database status
        """
        status = {}
        
        for genome, config in self.genome_configs.items():
            gtf_exists = config["gtf_file"].exists()
            db_exists = config["db_file"].exists()
            
            if db_exists:
                status[genome] = "ready"
            elif gtf_exists:
                status[genome] = "gtf_available"
            else:
                status[genome] = "not_downloaded"
        
        # Check for custom databases
        custom_dbs = list(self.data_dir.glob("custom_*.db"))
        for db_path in custom_dbs:
            genome_name = db_path.stem.replace("custom_", "")
            status[f"custom_{genome_name}"] = "ready"
        
        return status
    
    def get_database_info(self, genome: str) -> dict:
        """Get information about a specific genome database."""
        if genome in self.genome_configs:
            config = self.genome_configs[genome]
            return {
                "genome": genome,
                "gtf_url": config["gtf_url"],
                "gtf_file": str(config["gtf_file"]),
                "db_file": str(config["db_file"]),
                "gtf_exists": config["gtf_file"].exists(),
                "db_exists": config["db_file"].exists()
            }
        else:
            return {"error": f"Unknown genome: {genome}"}


# Global database manager instance
_db_manager = None

def get_db_manager() -> GeneDatabase:
    """Get global database manager instance."""
    global _db_manager
    if _db_manager is None:
        _db_manager = GeneDatabase()
    return _db_manager


def load_gene_database(genome: str = "hg38", gtf_file: Optional[str] = None) -> Optional[gffutils.FeatureDB]:
    """
    Convenience function to load gene database.
    
    Parameters:
    -----------
    genome : str, default "hg38"
        Genome assembly name
    gtf_file : str, optional
        Custom GTF file path (for genome="custom")
        
    Returns:
    --------
    gffutils.FeatureDB or None
        Loaded database
    """
    db_manager = get_db_manager()
    return db_manager.load_database(genome, gtf_file)