# Quick Start: Testing Pipeline in Dry-Run Mode

## What You Asked For

> "I would like to test the operation in the production env. first I would like to see what the pipeline would do. Just the backend. do not actually touch the data"

## What You Got

The pipeline now has a **dry-run mode** that shows exactly what would happen without modifying any data.

## How to Use It

### Option 1: Command Line (Simplest)

```bash
python pipeline.py --dry-run
```

That's it! The pipeline will:
- ✅ Show what operations would be performed
- ✅ Report what files would be processed
- ✅ Display all steps in detail
- ❌ NOT modify any files
- ❌ NOT write to the database

### Option 2: With Your Data Directory

```bash
python pipeline.py --dry-run --data-dir /path/to/your/data
```

### Option 3: Test Specific Phases Only

```bash
# Test only science target processing
python pipeline.py --dry-run --phase3-only

# Test only photometry setup
python pipeline.py --dry-run --phase4-only
```

## What You'll See

Output will look like this:

```
============================================================
DRY RUN MODE - No files or database will be modified
============================================================
Starting Wombat Pipeline
Data directory: ./data

[DRY RUN] Would connect to database and create tables if needed

Processing night: 20241115
  [DRY RUN] Would add Bias calibration: 10 images
  [DRY RUN] Would create master bias: mbias.fits
  [DRY RUN] Would add FlatV calibration: 15 images
  [DRY RUN] Would create master flat: mflatV.fits
  [DRY RUN] Would add target: V0450_Aql (filter: V, images: 50)

Processing science targets...
  [DRY RUN] Would create directory: ./data/20241115/Light/V0450_Aql/Reduced_images
  Calibrating: image001.fits
  [DRY RUN] Would save calibrated frame to: image001.fits
  Plate solving: image001.fits
  [DRY RUN] Would attempt plate solving with ASTAP for: image001.fits
  [DRY RUN] Would record WCS solution for image001.fits
...
```

Every line with `[DRY RUN]` is something that **would** happen but **won't** in dry-run mode.

## Files Created

Three new documentation files:

1. **DRY_RUN_README.md** - Complete documentation
2. **DRY_RUN_CHANGES.md** - Implementation details
3. **demo_dry_run.py** - Simple demo script

## Important Notes

- Your production data is **completely safe**
- The database is not touched (uses in-memory database instead)
- No FITS files are created or modified
- No directories are created
- No external programs (like ASTAP) are executed

## When You're Ready to Run For Real

Just remove the `--dry-run` flag:

```bash
# This WILL modify data
python pipeline.py --data-dir /path/to/your/data
```

## Need More Info?

See `DRY_RUN_README.md` for complete documentation.

---

**That's it!** You can now safely test the pipeline on your production environment without any risk to your data.
