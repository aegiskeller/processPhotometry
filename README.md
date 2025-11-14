# NightView - FWHM-based Differential Photometry

A Java/ImageJ application for performing differential photometry on FITS images with automatic FWHM-based aperture sizing.

## Features

- **FWHM-based aperture sizing** following AstroImageJ algorithm
  - Object aperture: 0.85 × FWHM  
  - Sky inner radius: 1.7 × FWHM
  - Sky outer radius: 3.4 × FWHM
- **Automatic star selection** via VSX and AAVSO VSP queries
- **Batch photometry** on multiple images
- **Interactive radial profile plots** for quality checking

## Quick Start

```bash
./build_and_run.sh
```

1. Select a target/night from the table
2. Image loads with detected stars overlay
3. Click stars to view radial profiles
4. Click "Do Phot All" to process all images in directory
5. Results appear in photometry table

## Project Structure

```
processPhotometry/
├── build_and_run.sh           # Launch script
├── target/classes/nightview/  # Compiled Java classes (working baseline)
├── build/jar/nightview.jar    # Application JAR
└── data/                      # FITS images by date
    └── YYYYMMDD/
        └── Light/
            └── TargetName/
```

## FWHM Algorithm

The application computes FWHM via radial profile analysis:

1. Calculate radial intensity profile from star center
2. Find peak intensity (center pixel)  
3. Determine background (mean of outer 1/3 of profile)
4. Calculate half-maximum = background + (peak - background) / 2
5. Find radius where intensity crosses half-maximum
6. FWHM = 2 × half-maximum radius

Apertures are then set automatically based on this FWHM measurement.

## Dependencies

- Java 11+
- ImageJ (ij.jar)
- Astronomy plugins (Astronomy_.jar)  
- JAMA matrix library (Jama.jar)

Libraries should be in `/Users/aegiskeller/Documents/lib/`

## Current Status

**Working:**
- ✅ FWHM-based aperture sizing
- ✅ VSX/VSP queries for star catalogs
- ✅ Differential photometry 
- ✅ Radial profile visualization
- ✅ Batch processing

**To Be Added:**
- Export photometry table to file
- Transformation to standard magnitudes
- AAVSO report generation

## Development

The application is built from pre-compiled class files in `target/classes`. To rebuild the JAR:

```bash
jar -cvf build/jar/nightview.jar -C target/classes .
```

Source code is decompiled and available in `src/` for reference, but cannot be recompiled due to decompiler artifacts.

## License

Open source for astronomical research.
