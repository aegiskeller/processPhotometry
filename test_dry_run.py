#!/usr/bin/env python3
"""
Test script to demonstrate dry-run mode for the pipeline.
"""

from pipeline import WombatPipeline
from pathlib import Path

def main():
    print("="*70)
    print("PIPELINE DRY-RUN TEST")
    print("="*70)
    print()
    print("This script will run the pipeline in DRY-RUN mode.")
    print("It will show what operations would be performed without:")
    print("  - Creating or modifying any files")
    print("  - Writing to the database")
    print("  - Moving or renaming files")
    print("  - Running plate solving")
    print()
    
    # Check if data directory exists
    data_dir = Path('./data')
    if not data_dir.exists():
        print(f"Warning: Data directory '{data_dir}' does not exist.")
        print("The dry-run will still show what would happen if it existed.")
        print()
    
    # Initialize pipeline in dry-run mode
    pipeline = WombatPipeline(
        data_dir='./data',
        db_path='wombatpipeline.db',
        dry_run=True  # This is the key parameter
    )
    
    try:
        # Run the full pipeline
        pipeline.run()
    finally:
        pipeline.close()
    
    print()
    print("="*70)
    print("DRY-RUN COMPLETE")
    print("="*70)
    print()
    print("Summary:")
    print("  - No files were created or modified")
    print("  - No database changes were made")
    print("  - All operations above were simulated only")
    print()

if __name__ == '__main__':
    main()
