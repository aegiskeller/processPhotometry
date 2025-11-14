#!/usr/bin/env python3
"""Test harness to examine FITS image pixel values and do aperture photometry."""

import numpy as np
from astropy.io import fits
import sys

def analyze_fits_image(filename):
    """Analyze a FITS image and report pixel statistics."""
    print(f"\n{'='*80}")
    print(f"Analyzing: {filename}")
    print(f"{'='*80}")
    
    with fits.open(filename) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        
        # Image info
        print(f"\nImage Information:")
        print(f"  BITPIX: {header.get('BITPIX', 'N/A')}")
        print(f"  Data type: {data.dtype}")
        print(f"  Shape: {data.shape}")
        
        # Pixel statistics
        print(f"\nPixel Value Statistics:")
        print(f"  Min:  {np.min(data):.2f}")
        print(f"  Max:  {np.max(data):.2f}")
        print(f"  Mean: {np.mean(data):.2f}")
        print(f"  Median: {np.median(data):.2f}")
        print(f"  Std:  {np.std(data):.2f}")
        
        # Percentiles
        print(f"\nPercentiles:")
        print(f"  1%:   {np.percentile(data, 1):.2f}")
        print(f"  25%:  {np.percentile(data, 25):.2f}")
        print(f"  50%:  {np.percentile(data, 50):.2f}")
        print(f"  75%:  {np.percentile(data, 75):.2f}")
        print(f"  99%:  {np.percentile(data, 99):.2f}")
        print(f"  99.9%: {np.percentile(data, 99.9):.2f}")
        
        # Find brightest pixel location
        max_loc = np.unravel_index(np.argmax(data), data.shape)
        print(f"\nBrightest pixel:")
        print(f"  Location (y,x): {max_loc}")
        print(f"  Value: {data[max_loc]:.2f}")
        
        # Sample around brightest pixel
        y, x = max_loc
        if y >= 5 and x >= 5 and y < data.shape[0]-5 and x < data.shape[1]-5:
            region = data[y-5:y+6, x-5:x+6]
            print(f"\n11x11 region around brightest pixel:")
            print(f"  Min:  {np.min(region):.2f}")
            print(f"  Max:  {np.max(region):.2f}")
            print(f"  Mean: {np.mean(region):.2f}")
        
        # Look for TX Eri coordinates if WCS is available
        if 'CRVAL1' in header and 'CD1_1' in header:
            print(f"\nWCS Information:")
            print(f"  CRVAL1 (RA):  {header['CRVAL1']:.6f}")
            print(f"  CRVAL2 (Dec): {header['CRVAL2']:.6f}")
            
            # TX Eri coordinates: RA=48.144958°, Dec=-19.532667°
            # Try to find TX Eri pixel location
            try:
                from astropy.wcs import WCS
                wcs = WCS(header)
                
                # TX Eri J2000 coordinates
                tx_eri_ra = 48.144958  # degrees
                tx_eri_dec = -19.532667  # degrees
                
                # Convert to pixel coordinates
                px_arr, py_arr = wcs.world_to_pixel_values(tx_eri_ra, tx_eri_dec)
                px, py = int(float(px_arr)), int(float(py_arr))
                
                print(f"\nTX Eri location:")
                print(f"  Catalog RA/Dec: {tx_eri_ra:.6f}, {tx_eri_dec:.6f}")
                print(f"  Pixel (x,y): ({px}, {py})")
                
                if 0 <= py < data.shape[0] and 0 <= px < data.shape[1]:
                    print(f"  Pixel value at TX Eri: {data[py, px]:.2f}")
                    
                    # Aperture photometry at TX Eri location
                    if py >= 15 and px >= 15 and py < data.shape[0]-15 and px < data.shape[1]-15:
                        print(f"\nAperture Photometry at TX Eri:")
                        
                        # Create circular aperture (radius 11 pixels)
                        aperture_radius = 11
                        y_grid, x_grid = np.ogrid[-aperture_radius:aperture_radius+1, 
                                                   -aperture_radius:aperture_radius+1]
                        aperture_mask = x_grid*x_grid + y_grid*y_grid <= aperture_radius*aperture_radius
                        
                        aperture_data = data[py-aperture_radius:py+aperture_radius+1, 
                                           px-aperture_radius:px+aperture_radius+1]
                        star_pixels = aperture_data[aperture_mask]
                        
                        print(f"  Aperture radius: {aperture_radius} pixels")
                        print(f"  Star pixels in aperture: {len(star_pixels)}")
                        print(f"  Min in aperture: {np.min(star_pixels):.2f}")
                        print(f"  Max in aperture: {np.max(star_pixels):.2f}")
                        print(f"  Mean in aperture: {np.mean(star_pixels):.2f}")
                        print(f"  Sum in aperture: {np.sum(star_pixels):.2f}")
                        
                        # Sky annulus (16-26 pixels)
                        sky_inner = 16
                        sky_outer = 26
                        y_grid_outer, x_grid_outer = np.ogrid[-sky_outer:sky_outer+1, 
                                                               -sky_outer:sky_outer+1]
                        y_grid_inner, x_grid_inner = np.ogrid[-sky_inner:sky_inner+1, 
                                                               -sky_inner:sky_inner+1]
                        sky_mask = ((x_grid_outer*x_grid_outer + y_grid_outer*y_grid_outer <= sky_outer*sky_outer) &
                                   (x_grid_inner*x_grid_inner + y_grid_inner*y_grid_inner > sky_inner*sky_inner))
                        
                        sky_data = data[py-sky_outer:py+sky_outer+1, 
                                      px-sky_outer:px+sky_outer+1]
                        sky_pixels = sky_data[sky_mask]
                        
                        print(f"\n  Sky annulus: {sky_inner}-{sky_outer} pixels")
                        print(f"  Sky pixels: {len(sky_pixels)}")
                        print(f"  Sky min: {np.min(sky_pixels):.2f}")
                        print(f"  Sky max: {np.max(sky_pixels):.2f}")
                        print(f"  Sky median: {np.median(sky_pixels):.2f}")
                        print(f"  Sky mean: {np.mean(sky_pixels):.2f}")
                        
                        # Mode approximation: 3*median - 2*mean
                        sky_mode = 3 * np.median(sky_pixels) - 2 * np.mean(sky_pixels)
                        print(f"  Sky mode (approx): {sky_mode:.2f}")
                        
                        # Net flux
                        star_sum = np.sum(star_pixels)
                        net_flux = star_sum - sky_mode * len(star_pixels)
                        print(f"\n  Star sum: {star_sum:.2f}")
                        print(f"  Sky contribution: {sky_mode * len(star_pixels):.2f}")
                        print(f"  Net flux: {net_flux:.2f}")
                        
                        if net_flux > 0:
                            inst_mag = -2.5 * np.log10(net_flux)
                            print(f"  Instrumental magnitude: {inst_mag:.3f}")
                        
            except ImportError:
                print("  (astropy.wcs not available for WCS analysis)")
            except Exception as e:
                print(f"  WCS analysis failed: {e}")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: test_image_values.py <fits_file> [<fits_file> ...]")
        print("\nTesting with default files...")
        files = [
            "data/20251109/Light/RRHya/2025-11-11_23-42-18_V_-10.00_73.00s_0062.fits",
            "data/20251109/Light/RRHya/Reduced_images/2025-11-11_23-42-18_V_-10.00_73.00s_0062_cal.fits",
            "data/20251111/Light/EETel/Reduced_images/2025-11-11_23-42-18_V_-10.00_73.00s_0062_cal.fits"
        ]
    else:
        files = sys.argv[1:]
    
    for filename in files:
        try:
            analyze_fits_image(filename)
        except FileNotFoundError:
            print(f"\nFile not found: {filename}")
        except Exception as e:
            print(f"\nError analyzing {filename}: {e}")
            import traceback
            traceback.print_exc()
    
    print(f"\n{'='*80}\n")
