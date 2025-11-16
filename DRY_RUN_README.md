# Dry-Run Mode for Production Testing

## Overview

The pipeline now supports a **dry-run mode** that allows you to test operations in your production environment without modifying any data. This is useful for:

- Validating pipeline logic before processing real data
- Testing with production data structure without risk
- Debugging pipeline behavior
- Documenting what the pipeline would do

## What Dry-Run Mode Does

When running in dry-run mode, the pipeline will:

✅ **Show what would happen:**
- Print all operations that would be performed
- Display file paths that would be created or modified
- Show database operations that would be executed
- Report what data would be processed

❌ **NOT actually do:**
- Create or modify any FITS files
- Write to the database (uses in-memory database instead)
- Move or rename files
- Create directories
- Run plate solving with ASTAP
- Write output files (photometry tables, reports, etc.)

## How to Use Dry-Run Mode

### Method 1: Command Line (Recommended)

Run the pipeline with the `--dry-run` flag:

```bash
python pipeline.py --dry-run
```

Additional options:

```bash
# Dry-run with custom data directory
python pipeline.py --dry-run --data-dir /path/to/data

# Dry-run only for science target processing (Phase 3)
python pipeline.py --dry-run --phase3-only

# Dry-run only for photometry setup (Phase 4)
python pipeline.py --dry-run --phase4-only
```

### Method 2: Test Script

Use the provided test script:

```bash
python test_dry_run.py
```

### Method 3: In Code

Import and use programmatically:

```python
from pipeline import WombatPipeline

# Initialize with dry_run=True
pipeline = WombatPipeline(
    data_dir='./data',
    db_path='wombatpipeline.db',
    dry_run=True
)

try:
    pipeline.run()
finally:
    pipeline.close()
```

## Understanding Dry-Run Output

Look for these markers in the output:

- `[DRY RUN]` - Indicates a simulated operation
- `Would...` - Describes what would happen in a real run

Example output:
```
[DRY RUN] Would create directory: /path/to/Light/target/Reduced_images
[DRY RUN] Would add Bias calibration: 10 images
[DRY RUN] Would create master bias: mbias.fits
[DRY RUN] Would save calibrated frame to: image001.fits
[DRY RUN] Would attempt plate solving with ASTAP for: image001.fits
[DRY RUN] Would record WCS solution for image001.fits
[DRY RUN] Would update target table with proc_bias=20241115, proc_flat=20241115
```

## Complete Command Reference

```bash
# Full help
python pipeline.py --help

# Basic dry-run
python pipeline.py --dry-run

# Dry-run with specific data directory
python pipeline.py --dry-run --data-dir /mnt/observatory/data

# Dry-run with custom database path (not written in dry-run mode)
python pipeline.py --dry-run --db-path custom.db

# Phase-specific dry-runs
python pipeline.py --dry-run --phase3-only  # Only science targets
python pipeline.py --dry-run --phase4-only  # Only photometry setup
```

## What Gets Checked

Even in dry-run mode, the pipeline will:

1. **Validate directory structure** - Check if expected directories exist
2. **Read FITS headers** - Extract filter information, coordinates, etc.
3. **Query databases** - VSX lookups, AAVSO VSP queries (read-only)
4. **Count files** - Report how many files would be processed
5. **Check calibrations** - Find closest bias/flat frames
6. **Simulate processing** - Show step-by-step what would happen

## Safety Features

- Uses in-memory SQLite database (`:memory:`) instead of the real database
- All file operations are logged but not executed
- Original data remains completely untouched
- No external commands are executed (e.g., ASTAP)

## Best Practices

1. **Always dry-run first** in production before processing new data
2. **Review the output** to ensure expected behavior
3. **Check file counts** match expectations
4. **Verify directory paths** are correct
5. **Test each phase** separately if needed

## Troubleshooting

### "Data directory does not exist"
- The dry-run will still run but won't find any data to process
- Check your `--data-dir` path

### "No Light directory found"
- Normal if testing with empty data structure
- Indicates which nights would be skipped

### Plate solving shows "Skipped in dry-run mode"
- Expected behavior - ASTAP is not called in dry-run mode
- Real run would attempt actual plate solving

## Example Workflow

```bash
# 1. Test on production data without changes
python pipeline.py --dry-run --data-dir /data/observatory/2024-11-15

# 2. Review output, verify expected behavior

# 3. If everything looks good, run for real (remove --dry-run)
python pipeline.py --data-dir /data/observatory/2024-11-15
```

## Technical Details

### Modified Operations

All these operations are simulated in dry-run mode:

- `mkdir()` - Directory creation
- `rename()` - File moves/renames  
- `writeto()` - FITS file writing
- `open(file, 'w')` - Text file writing
- `cursor.execute(INSERT/UPDATE)` - Database modifications
- `subprocess.run(astap)` - External program execution

### Database Behavior

In dry-run mode:
- Creates in-memory database (`sqlite3.connect(':memory:')`)
- All schema creation happens normally
- INSERT/UPDATE/DELETE operations are logged but use memory
- Database is discarded when pipeline exits

This allows the pipeline to maintain state during the dry-run for testing logic without persisting any data.
