# Pipeline Dry-Run Mode Implementation Summary

## Overview

The pipeline backend (`pipeline.py`) has been enhanced with a comprehensive dry-run mode that allows you to test operations in your production environment without modifying any data.

## Files Modified

### 1. `pipeline.py` (Main Pipeline)

**Key Changes:**

- **Constructor (`__init__`)**: Added `dry_run` parameter (default: `False`)
  - Displays banner when dry-run mode is enabled
  - Stores dry_run flag for use throughout the class

- **Database Operations**: Modified to check `self.dry_run` before executing:
  - `setup_database()` - Uses in-memory database instead of file
  - All `INSERT` statements - Logged but not executed
  - All `UPDATE` statements - Logged but not executed
  - `update_target_status()` - Early return with log message

- **File Operations**: Protected against modifications:
  - `organize_flats_by_filter()` - Directory creation and file moves logged only
  - `create_master_bias()` - FITS writing skipped
  - `create_master_flat()` - FITS writing skipped
  - `calibrate_frame()` - Output file writing skipped
  - `_write_photometry_table()` - File writing skipped
  - `create_aavso_report()` - File writing skipped
  - `process_phase6_transformation()` - Transformed table writing skipped
  - `process_science_targets()` - Reduced_images directory creation logged only

- **External Program Execution**:
  - `plate_solve_astap()` - Returns mock failure result without calling ASTAP

- **Command-Line Interface**: Added argument parsing:
  - `--dry-run` - Enable dry-run mode
  - `--data-dir` - Specify data directory
  - `--db-path` - Specify database path
  - `--phase3-only` - Process only science targets
  - `--phase4-only` - Process only photometry setup

## Files Created

### 2. `DRY_RUN_README.md`

Complete documentation including:
- What dry-run mode does and doesn't do
- How to use it (3 different methods)
- Command reference with examples
- Understanding output
- Best practices
- Troubleshooting guide
- Technical implementation details

### 3. `test_dry_run.py`

Standalone test script that:
- Demonstrates dry-run functionality
- Provides clear user feedback
- Runs full pipeline in dry-run mode
- Shows summary of what was (not) done

### 4. `demo_dry_run.py`

Educational demonstration script that:
- Explains the changes made
- Shows usage examples
- Lists all protected operations
- Provides quick reference

## How Dry-Run Mode Works

### Initialization
```python
pipeline = WombatPipeline(dry_run=True)
```

### Database Protection
- Creates in-memory database: `sqlite3.connect(':memory:')`
- All writes logged but executed only in memory
- Real database remains untouched

### File Operation Protection
Every destructive operation checks:
```python
if self.dry_run:
    print(f"[DRY RUN] Would {operation_description}")
    return
# actual operation here
```

### Output Markers
All dry-run operations are clearly marked:
- `[DRY RUN]` prefix for simulated operations
- `Would...` phrasing to indicate future tense

## Usage Examples

### Basic Dry-Run
```bash
python pipeline.py --dry-run
```

### With Custom Data Directory
```bash
python pipeline.py --dry-run --data-dir /mnt/observatory/data/20241115
```

### Phase-Specific Testing
```bash
python pipeline.py --dry-run --phase3-only  # Science targets only
python pipeline.py --dry-run --phase4-only  # Photometry only
```

### Programmatic Usage
```python
from pipeline import WombatPipeline

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

## Protected Operations

### Database Writes
- INSERT INTO cal (calibration data)
- INSERT INTO target (target data)
- INSERT OR REPLACE INTO wcs (plate solving results)
- UPDATE target (status updates)

### File Creation
- Master bias FITS files (mbias.fits)
- Master flat FITS files (mflatX.fits)
- Calibrated science frames (Reduced_images/*.fits)
- Photometry tables (*.tbl)
- Transformed photometry tables (*_transformed.tbl)
- AAVSO reports (AAVSO_AID_Report_*.txt)

### File Modifications
- FITS header updates (plate solving)
- File moves (flat organization)

### Directory Creation
- Filter subdirectories (FlatV, FlatB, etc.)
- Reduced_images directories

### External Programs
- ASTAP plate solving

## Testing Performed

✅ Syntax validation: `python -m py_compile pipeline.py`
✅ Demo script execution: `python demo_dry_run.py`
✅ Code review of all file operations
✅ Code review of all database operations

## Recommendations

1. **Always test with dry-run first** on production data
2. **Review the output** to verify expected behavior
3. **Install dependencies** to run actual tests:
   ```bash
   pip install -r requirements.txt  # if available
   ```
4. **Run on real data structure** to see actual processing flow
5. **Compare dry-run output** with expected results before real run

## Next Steps

To use this in production:

1. Install all required dependencies (astropy, photutils, etc.)
2. Run dry-run on a test dataset:
   ```bash
   python pipeline.py --dry-run --data-dir /path/to/test/data
   ```
3. Review output carefully
4. When satisfied, remove `--dry-run` flag to process for real

## Safety Features

- ✅ In-memory database prevents real DB writes
- ✅ All file operations gated by dry_run check
- ✅ External programs not executed
- ✅ Clear output markers ([DRY RUN])
- ✅ Summary at end confirms no changes made
- ✅ Original data completely untouched

## Benefits

1. **Risk-free testing** - Test with real data without consequences
2. **Validation** - Verify pipeline logic before running
3. **Documentation** - Generate processing plan for review
4. **Debugging** - Understand what pipeline would do
5. **Education** - Learn how pipeline processes data
