# Quick Start - Differential Photometry Workflow

## Complete Workflow

### Step 1: Aperture Photometry

```bash
./build_and_run.sh
```

1. Select target/night from the table  
2. Stars are automatically detected and labeled
3. Click "Do Phot All" to process all images
4. **Manual step**: Right-click the table → "Save As..." → `<TargetName>-<Date>.tbl`

### Step 2: Transform to Standard Magnitudes  

```bash
python3 analyze_photometry.py <TargetName>-<Date>.tbl --fetch-aavso --transform
```

This will:
- Query AAVSO for comparison star standard magnitudes
- Match your stars to catalog  
- Compute transformation coefficients
- Apply to target and check star

### Step 3: View Results (Optional GUI)

```bash
python3 photometry_gui.py
```

Load your .tbl file, view time series plots, and export AAVSO reports.

## Example

```bash
# Run photometry
./build_and_run.sh
# (Select RRLyr, click "Do Phot All", save table as RRLyr-2025-11-09.tbl)

# Transform
python3 analyze_photometry.py RRLyr-2025-11-09.tbl --fetch-aavso --transform

# GUI analysis
python3 photometry_gui.py
```

## What's Happening

**Aperture Photometry:**
- FWHM calculated for each star via radial profile
- Apertures set to 0.85×, 1.7×, 3.4× FWHM
- Instrumental magnitudes computed
- Results saved in tab-delimited format

**Transformation:**
- Comparison stars matched to AAVSO catalog (2" tolerance)
- Linear fit: V_std = V_inst + ZP + k×(V-I)
- Applied to target (T1) and check (T2)
- Quality metrics computed (RMS, scatter)

## Quality Checks

- **RMS < 0.02 mag**: Excellent transformation
- **Check star scatter < 0.03 mag**: Stable conditions
- **Time series plot**: Target should vary, check should be flat

## Files

- Input: `data/YYYYMMDD/Light/TargetName/*.fits`
- Output: `<TargetName>-<Date>.tbl` (photometry table)
- Final: AAVSO Extended Format report

---

For detailed information, see README.md
