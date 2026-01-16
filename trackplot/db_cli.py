"""
Database management CLI for trackplot package.
"""

import argparse
import sys
from .database import get_db_manager


def main():
    """Database management CLI."""
    parser = argparse.ArgumentParser(
        description="Manage trackplot gene databases",
        prog="trackplot-db"
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # List command
    list_parser = subparsers.add_parser('list', help='List available databases')
    
    # Download command  
    download_parser = subparsers.add_parser('download', help='Download genome database')
    download_parser.add_argument('genome', choices=['hg38', 'hg19', 'mm10', 'mm39'],
                               help='Genome assembly to download')
    download_parser.add_argument('--force', action='store_true',
                               help='Force re-download even if exists')
    
    # Info command
    info_parser = subparsers.add_parser('info', help='Show database information')
    info_parser.add_argument('genome', help='Genome assembly name')
    
    # Create command
    create_parser = subparsers.add_parser('create', help='Create database from custom GTF')
    create_parser.add_argument('gtf_file', help='Path to GTF file')
    create_parser.add_argument('--name', help='Custom database name')
    create_parser.add_argument('--force', action='store_true',
                              help='Force re-creation even if exists')
    
    args = parser.parse_args()
    
    if args.command is None:
        parser.print_help()
        return
    
    db_manager = get_db_manager()
    
    if args.command == 'list':
        print("Available genome databases:")
        print("-" * 40)
        status = db_manager.list_available_databases()
        for genome, state in status.items():
            status_icon = {
                'ready': '✅',
                'gtf_available': '📄', 
                'not_downloaded': '❌'
            }.get(state, '❓')
            print(f"{status_icon} {genome:15} {state}")
    
    elif args.command == 'download':
        try:
            db_file = db_manager.create_database(args.genome, force=args.force)
            print(f"✅ Database ready: {db_file}")
        except Exception as e:
            print(f"❌ Error: {e}")
            sys.exit(1)
    
    elif args.command == 'info':
        info = db_manager.get_database_info(args.genome)
        if 'error' in info:
            print(f"❌ {info['error']}")
        else:
            print(f"Genome: {info['genome']}")
            print(f"GTF URL: {info['gtf_url']}")
            print(f"GTF file: {info['gtf_file']} ({'exists' if info['gtf_exists'] else 'missing'})")
            print(f"Database: {info['db_file']} ({'exists' if info['db_exists'] else 'missing'})")
    
    elif args.command == 'create':
        try:
            name = args.name or f"custom_{Path(args.gtf_file).stem}"
            db_file = db_manager.create_database("custom", args.gtf_file, force=args.force)
            print(f"✅ Custom database created: {db_file}")
        except Exception as e:
            print(f"❌ Error: {e}")
            sys.exit(1)


if __name__ == "__main__":
    main()