#!/usr/bin/env python3
"""
Automated astronomical data processing pipeline backend.
Scans data directory, organizes calibration files by filter, and populates SQLite database.
"""

import os
import sqlite3
import glob
from pathlib import Path
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u
import numpy as np
import subprocess
import re
from datetime import datetime
import requests
from urllib.parse import urlencode


class WombatPipeline:
    def __init__(self, data_dir='./data', db_path='wombatpipeline.db', astap_path='/Applications/ASTAP.app/Contents/MacOS/astap'):
        """
        Initialize the pipeline.
        
        Args:
            data_dir: Root directory containing night subdirectories
            db_path: Path to SQLite database
            astap_path: Path to ASTAP executable
        """
        self.data_dir = Path(data_dir)
        self.db_path = db_path
        self.astap_path = astap_path
        self.conn = None
        
    def setup_database(self):
        """Create database and tables if they don't exist."""
        self.conn = sqlite3.connect(self.db_path)
        cursor = self.conn.cursor()
        
        # Create target table
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS target (
                date TEXT NOT NULL,
                target_name TEXT NOT NULL,
                filter TEXT,
                num_images INTEGER DEFAULT 0,
                proc_bias INTEGER DEFAULT 0,
                proc_flat INTEGER DEFAULT 0,
                apphot INTEGER DEFAULT 0,
                submitted INTEGER DEFAULT 0,
                PRIMARY KEY (date, target_name, filter)
            )
        ''')
        
        # Create calibration table
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS cal (
                date TEXT NOT NULL,
                cal_type TEXT NOT NULL,
                filter TEXT,
                num_images INTEGER,
                PRIMARY KEY (date, cal_type, filter)
            )
        ''')
        
        # Create WCS table for plate solving results
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS wcs (
                date TEXT NOT NULL,
                target_name TEXT NOT NULL,
                filter TEXT,
                filename TEXT NOT NULL,
                success INTEGER DEFAULT 0,
                ra_deg REAL,
                dec_deg REAL,
                solve_time REAL,
                num_stars INTEGER,
                faintest_mag REAL,
                fov_deg REAL,
                error_message TEXT,
                timestamp TEXT,
                PRIMARY KEY (date, target_name, filter, filename)
            )
        ''')
        
        self.conn.commit()
        print(f"Database setup complete: {self.db_path}")
        
    def get_filter_from_fits(self, fits_path):
        """
        Extract filter information from FITS header.
        
        Args:
            fits_path: Path to FITS file
            
        Returns:
            Filter name or None if not found
        """
        try:
            with fits.open(fits_path) as hdul:
                header = hdul[0].header
                # Try common filter header keywords
                for key in ['FILTER', 'FILTER1', 'FILTNAM', 'FILT']:
                    if key in header:
                        return header[key].strip()
        except Exception as e:
            print(f"Warning: Could not read filter from {fits_path}: {e}")
        return None
        
    def organize_flats_by_filter(self, flat_dir):
        """
        Organize flat frames into subdirectories by filter.
        
        Args:
            flat_dir: Path to Flat directory
            
        Returns:
            Dictionary mapping filters to list of FITS files
        """
        filter_files = {}
        fits_files = list(flat_dir.glob('*.fits'))
        
        for fits_file in fits_files:
            filter_name = self.get_filter_from_fits(fits_file)
            if filter_name:
                if filter_name not in filter_files:
                    filter_files[filter_name] = []
                filter_files[filter_name].append(fits_file)
        
        # Create subdirectories and move files
        for filter_name, files in filter_files.items():
            filter_subdir = flat_dir / f"Flat{filter_name}"
            filter_subdir.mkdir(exist_ok=True)
            
            for fits_file in files:
                dest = filter_subdir / fits_file.name
                if not dest.exists():
                    fits_file.rename(dest)
                    
        return filter_files
        
    def count_images(self, directory):
        """Count FITS images in a directory."""
        if not directory.exists():
            return 0
        return len(list(directory.glob('*.fits')))
        
    def create_master_bias(self, bias_dir):
        """
        Create master bias from individual bias frames using median combination with sigma clipping.
        
        Args:
            bias_dir: Path to Bias directory
            
        Returns:
            True if master bias was created, False otherwise
        """
        mbias_path = bias_dir / 'mbias.fits'
        
        # Skip if master bias already exists
        if mbias_path.exists():
            print(f"    Master bias already exists - skipping")
            return False
            
        # Get all bias frames
        bias_files = sorted(bias_dir.glob('*.fits'))
        if len(bias_files) < 3:
            print(f"    Not enough bias frames ({len(bias_files)} < 3) - skipping master bias creation")
            return False
            
        print(f"    Creating master bias from {len(bias_files)} frames...")
        
        # Load all bias frames into a 3D array
        bias_data = []
        for bias_file in bias_files:
            with fits.open(bias_file) as hdul:
                bias_data.append(hdul[0].data.astype(np.float32))
        
        bias_cube = np.array(bias_data)
        
        # Median combine with sigma clipping
        # Clip outliers at 3 sigma, then take median
        median_bias = np.median(bias_cube, axis=0)
        std_bias = np.std(bias_cube, axis=0)
        
        # Mask outliers
        mask = np.abs(bias_cube - median_bias) > 3 * std_bias
        masked_cube = np.ma.array(bias_cube, mask=mask)
        
        # Final median
        master_bias = np.ma.median(masked_cube, axis=0).filled(median_bias)
        
        # Save master bias
        with fits.open(bias_files[0]) as hdul:
            hdu = fits.PrimaryHDU(master_bias.astype(np.float32), header=hdul[0].header)
            hdu.header['HISTORY'] = f'Master bias created from {len(bias_files)} frames'
            hdu.header['HISTORY'] = 'Median combined with 3-sigma clipping'
            hdu.writeto(mbias_path, overwrite=True)
            
        print(f"    Created master bias: {mbias_path.name}")
        return True
        
    def create_master_flat(self, flat_dir, filter_name):
        """
        Create master flat from individual flat frames using normalized median combination with sigma clipping.
        
        Args:
            flat_dir: Path to FlatX directory (e.g., FlatV, FlatB)
            filter_name: Filter name (e.g., 'V', 'B')
            
        Returns:
            True if master flat was created, False otherwise
        """
        mflat_path = flat_dir / f'mflat{filter_name}.fits'
        
        # Skip if master flat already exists
        if mflat_path.exists():
            print(f"    Master flat already exists - skipping")
            return False
            
        # Get all flat frames
        flat_files = sorted(flat_dir.glob('*.fits'))
        # Exclude any existing master flat
        flat_files = [f for f in flat_files if not f.name.startswith('mflat')]
        
        if len(flat_files) < 3:
            print(f"    Not enough flat frames ({len(flat_files)} < 3) - skipping master flat creation")
            return False
            
        print(f"    Creating master flat{filter_name} from {len(flat_files)} frames...")
        
        # Load and normalize all flat frames
        normalized_flats = []
        for flat_file in flat_files:
            with fits.open(flat_file) as hdul:
                flat_data = hdul[0].data.astype(np.float32)
                # Normalize by median
                median_val = np.median(flat_data)
                if median_val > 0:
                    normalized_flats.append(flat_data / median_val)
                else:
                    print(f"    Warning: {flat_file.name} has median <= 0, skipping")
        
        if len(normalized_flats) < 3:
            print(f"    Not enough valid flat frames after normalization - skipping")
            return False
            
        flat_cube = np.array(normalized_flats)
        
        # Median combine with sigma clipping
        median_flat = np.median(flat_cube, axis=0)
        std_flat = np.std(flat_cube, axis=0)
        
        # Mask outliers at 3 sigma
        mask = np.abs(flat_cube - median_flat) > 3 * std_flat
        masked_cube = np.ma.array(flat_cube, mask=mask)
        
        # Final median
        master_flat = np.ma.median(masked_cube, axis=0).filled(median_flat)
        
        # Save master flat
        with fits.open(flat_files[0]) as hdul:
            hdu = fits.PrimaryHDU(master_flat.astype(np.float32), header=hdul[0].header)
            hdu.header['HISTORY'] = f'Master flat created from {len(normalized_flats)} frames'
            hdu.header['HISTORY'] = 'Normalized and median combined with 3-sigma clipping'
            hdu.writeto(mflat_path, overwrite=True)
            
        print(f"    Created master flat: {mflat_path.name}")
        return True
        
    def process_night(self, night_dir):
        """
        Process a single night directory.
        
        Args:
            night_dir: Path to night directory
        """
        night_name = night_dir.name
        print(f"\nProcessing night: {night_name}")
        
        bias_dir = night_dir / 'Bias'
        flat_dir = night_dir / 'Flat'
        light_dir = night_dir / 'Light'
        
        # Check if Light directory exists - if not, ignore this night
        if not light_dir.exists():
            print(f"  No Light directory found - skipping night {night_name}")
            return
            
        cursor = self.conn.cursor()
        
        # Process Bias calibrations
        if bias_dir.exists():
            num_bias = self.count_images(bias_dir)
            if num_bias > 0:
                cursor.execute('''
                    INSERT OR REPLACE INTO cal (date, cal_type, filter, num_images)
                    VALUES (?, ?, ?, ?)
                ''', (night_name, 'Bias', None, num_bias))
                print(f"  Added Bias calibration: {num_bias} images")
                
                # Create master bias if needed
                self.create_master_bias(bias_dir)
        
        # Process Flat calibrations
        if flat_dir.exists():
            # First organize flats by filter if needed
            fits_in_flat = list(flat_dir.glob('*.fits'))
            if fits_in_flat:
                print(f"  Organizing {len(fits_in_flat)} flat frames by filter...")
                filter_files = self.organize_flats_by_filter(flat_dir)
            
            # Now count flats in each filter subdirectory
            for filter_subdir in flat_dir.glob('Flat*'):
                if filter_subdir.is_dir():
                    filter_name = filter_subdir.name.replace('Flat', '')
                    num_flats = self.count_images(filter_subdir)
                    if num_flats > 0:
                        cursor.execute('''
                            INSERT OR REPLACE INTO cal (date, cal_type, filter, num_images)
                            VALUES (?, ?, ?, ?)
                        ''', (night_name, 'Flat', filter_name, num_flats))
                        print(f"  Added Flat{filter_name} calibration: {num_flats} images")
                        
                        # Create master flat if needed
                        self.create_master_flat(filter_subdir, filter_name)
        
        # Process Light targets
        for target_dir in light_dir.iterdir():
            if target_dir.is_dir():
                target_name = target_dir.name
                
                # Get all FITS files and group by filter
                fits_files = list(target_dir.glob('*.fits'))
                filter_counts = {}
                
                for fits_file in fits_files:
                    filter_name = self.get_filter_from_fits(fits_file)
                    if filter_name not in filter_counts:
                        filter_counts[filter_name] = 0
                    filter_counts[filter_name] += 1
                
                # Insert a row for each filter
                for filter_name, num_images in filter_counts.items():
                    try:
                        cursor.execute('''
                            INSERT INTO target (date, target_name, filter, num_images, proc_bias, proc_flat, apphot, submitted)
                            VALUES (?, ?, ?, ?, 0, 0, 0, 0)
                        ''', (night_name, target_name, filter_name, num_images))
                        print(f"  Added target: {target_name} (filter: {filter_name}, images: {num_images})")
                    except sqlite3.IntegrityError:
                        print(f"  Target {target_name} with filter {filter_name} already exists for {night_name} - skipping")
        
        self.conn.commit()
        
    def find_closest_calibration(self, target_date, cal_type, filter_name=None):
        """
        Find the closest calibration file (forward or backward in time).
        
        Args:
            target_date: Date string of the target night
            cal_type: 'Bias' or 'Flat'
            filter_name: Filter name (required for Flat, None for Bias)
            
        Returns:
            Tuple of (cal_date, cal_path) or (None, None) if not found
        """
        cursor = self.conn.cursor()
        
        if cal_type == 'Bias':
            cursor.execute('''
                SELECT date FROM cal 
                WHERE cal_type = 'Bias' 
                ORDER BY ABS(JULIANDAY(date) - JULIANDAY(?))
                LIMIT 1
            ''', (target_date,))
        else:  # Flat
            cursor.execute('''
                SELECT date FROM cal 
                WHERE cal_type = 'Flat' AND filter = ?
                ORDER BY ABS(JULIANDAY(date) - JULIANDAY(?))
                LIMIT 1
            ''', (filter_name, target_date))
        
        result = cursor.fetchone()
        if not result:
            return None, None
            
        cal_date = result[0]
        
        # Build path to calibration file
        if cal_type == 'Bias':
            cal_path = self.data_dir / cal_date / 'Bias' / 'mbias.fits'
        else:
            cal_path = self.data_dir / cal_date / 'Flat' / f'Flat{filter_name}' / f'mflat{filter_name}.fits'
        
        if cal_path.exists():
            return cal_date, cal_path
        return None, None
        
    def calibrate_frame(self, science_path, bias_path=None, flat_path=None, output_path=None):
        """
        Calibrate a science frame with bias and flat.
        
        Args:
            science_path: Path to science FITS file
            bias_path: Path to master bias (or None)
            flat_path: Path to master flat (or None)
            output_path: Path for calibrated output
            
        Returns:
            True if successful, False otherwise
        """
        try:
            with fits.open(science_path) as hdul:
                science_data = hdul[0].data.astype(np.float32)
                header = hdul[0].header.copy()
                
                # Bias subtraction
                if bias_path and bias_path.exists():
                    with fits.open(bias_path) as bias_hdul:
                        bias_data = bias_hdul[0].data.astype(np.float32)
                        science_data = science_data - bias_data
                        header['HISTORY'] = f'Bias subtracted: {bias_path.name}'
                
                # Flat division
                if flat_path and flat_path.exists():
                    with fits.open(flat_path) as flat_hdul:
                        flat_data = flat_hdul[0].data.astype(np.float32)
                        # Avoid division by zero
                        flat_data[flat_data == 0] = 1.0
                        science_data = science_data / flat_data
                        header['HISTORY'] = f'Flat corrected: {flat_path.name}'
                
                # Save calibrated frame
                if output_path:
                    hdu = fits.PrimaryHDU(science_data, header=header)
                    hdu.writeto(output_path, overwrite=True)
                    
            return True
        except Exception as e:
            print(f"    Error calibrating {science_path.name}: {e}")
            return False
            
    def plate_solve_astap(self, fits_file, fov_attempts=[10, 30]):
        """
        Plate solve a FITS file using ASTAP.
        
        Args:
            fits_file: Path to FITS file
            fov_attempts: List of FOV values to try in degrees
            
        Returns:
            Dictionary with solution info
        """
        for fov in fov_attempts:
            cmd = [
                self.astap_path,
                '-f', str(fits_file),
                '-r', str(fov),
                '-update'  # Update FITS header with WCS solution
            ]
            
            try:
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=60
                )
                
                output = result.stdout + result.stderr
                
                if result.returncode == 0 and "Solution found:" in output:
                    # Parse solution
                    solution = {
                        'success': True,
                        'fov_deg': fov,
                        'timestamp': datetime.now().isoformat()
                    }
                    
                    # Extract RA/Dec
                    ra_match = re.search(r'Solution found:\s+(\d+):\s*(\d+)\s+([\d.]+)', output)
                    dec_match = re.search(r'(-?\d+)°\s*(\d+)\s+([\d.]+)', output)
                    
                    if ra_match:
                        ra_h, ra_m, ra_s = ra_match.groups()
                        solution['ra_deg'] = float(ra_h) * 15 + float(ra_m) / 4 + float(ra_s) / 240
                    
                    if dec_match:
                        dec_d, dec_m, dec_s = dec_match.groups()
                        dec_deg = abs(float(dec_d)) + float(dec_m) / 60 + float(dec_s) / 3600
                        if dec_d.startswith('-'):
                            dec_deg *= -1
                        solution['dec_deg'] = dec_deg
                    
                    # Extract solve time
                    time_match = re.search(r'Solved in ([\d.]+) sec', output)
                    if time_match:
                        solution['solve_time'] = float(time_match.group(1))
                    
                    # Extract number of stars
                    stars_match = re.search(r'(\d+) stars', output)
                    if stars_match:
                        solution['num_stars'] = int(stars_match.group(1))
                    
                    # Extract faintest magnitude
                    mag_match = re.search(r'magnitude: ([\d.]+)', output)
                    if mag_match:
                        solution['faintest_mag'] = float(mag_match.group(1))
                    
                    return solution
                    
            except subprocess.TimeoutExpired:
                print(f"    ASTAP timeout with FOV={fov}°")
            except Exception as e:
                print(f"    ASTAP error with FOV={fov}°: {e}")
        
        # All attempts failed
        return {
            'success': False,
            'error_message': f'Failed to solve with FOV attempts: {fov_attempts}',
            'timestamp': datetime.now().isoformat()
        }
        
    def query_vsx(self, target_name):
        """
        Query AAVSO VSX for target information.
        
        Args:
            target_name: Name of the variable star
            
        Returns:
            Dictionary with target info or None if not found
        """
        try:
            # VSX API endpoint
            url = "https://www.aavso.org/vsx/index.php"
            params = {
                'view': 'api.object',
                'format': 'json',
                'ident': target_name
            }
            
            response = requests.get(url, params=params, timeout=30)
            
            if response.status_code == 200:
                data = response.json()
                
                if 'VSXObject' in data and data['VSXObject']:
                    obj = data['VSXObject']
                    
                    # Parse RA/Dec
                    ra_str = obj.get('RA2000')
                    dec_str = obj.get('Declination2000')
                    
                    if ra_str and dec_str:
                        # Convert to decimal degrees
                        coord = SkyCoord(ra_str, dec_str, unit=(u.hourangle, u.deg))
                        
                        return {
                            'name': obj.get('Name'),
                            'auid': obj.get('AUID'),
                            'ra_deg': coord.ra.deg,
                            'dec_deg': coord.dec.deg,
                            'ra_hms': ra_str,
                            'dec_dms': dec_str,
                            'var_type': obj.get('Type'),
                            'period': obj.get('Period'),
                            'mag_max': obj.get('MaxMag'),
                            'mag_min': obj.get('MinMag'),
                            'epoch': obj.get('Epoch')
                        }
            
            return None
            
        except Exception as e:
            print(f"    Error querying VSX for {target_name}: {e}")
            return None
            
    def calculate_time_and_airmass(self, fits_file, target_ra, target_dec, observatory_location=None):
        """
        Calculate observation time (JD, HJD) and airmass from FITS header.
        
        Args:
            fits_file: Path to FITS file
            target_ra: Target RA in degrees
            target_dec: Target Dec in degrees
            observatory_location: EarthLocation object (default: use header or estimate)
            
        Returns:
            Dictionary with time and airmass info
        """
        try:
            with fits.open(fits_file) as hdul:
                header = hdul[0].header
                
                # Get observation time from header
                if 'DATE-OBS' in header:
                    obs_time_str = header['DATE-OBS']
                    # Handle various formats
                    if 'T' not in obs_time_str and 'TIME-OBS' in header:
                        obs_time_str = obs_time_str + 'T' + header['TIME-OBS']
                elif 'JD' in header:
                    obs_time = Time(header['JD'], format='jd')
                else:
                    print(f"    Warning: No time information in FITS header")
                    return None
                
                if 'obs_time_str' in locals():
                    obs_time = Time(obs_time_str, format='isot', scale='utc')
                
                # Get or estimate observatory location
                if observatory_location is None:
                    # Try to get from header
                    if all(k in header for k in ['SITELAT', 'SITELONG', 'SITEELEV']):
                        observatory_location = EarthLocation(
                            lat=header['SITELAT'] * u.deg,
                            lon=header['SITELONG'] * u.deg,
                            height=header['SITEELEV'] * u.m
                        )
                    else:
                        # Use a default location (needs to be configured)
                        print(f"    Warning: No site coordinates in header, using default")
                        # Default to Greenwich as placeholder
                        observatory_location = EarthLocation(lat=51.4778*u.deg, lon=0*u.deg, height=0*u.m)
                
                # Calculate JD and HJD
                jd = obs_time.jd
                
                # Heliocentric correction
                target_coord = SkyCoord(ra=target_ra*u.deg, dec=target_dec*u.deg, frame='icrs')
                ltt_helio = obs_time.light_travel_time(target_coord, kind='heliocentric', location=observatory_location)
                hjd = (obs_time + ltt_helio).jd
                
                # Calculate airmass
                altaz_frame = AltAz(obstime=obs_time, location=observatory_location)
                target_altaz = target_coord.transform_to(altaz_frame)
                
                # Airmass using simple sec(z) approximation
                zenith_angle = 90 - target_altaz.alt.deg
                if zenith_angle < 90:
                    airmass = 1.0 / np.cos(np.radians(zenith_angle))
                else:
                    airmass = 99.0  # Below horizon
                
                return {
                    'jd': jd,
                    'hjd': hjd,
                    'airmass': airmass,
                    'altitude': target_altaz.alt.deg,
                    'azimuth': target_altaz.az.deg,
                    'obs_time': obs_time.isot
                }
                
        except Exception as e:
            print(f"    Error calculating time/airmass: {e}")
            return None
            
    def get_fov_from_wcs(self, fits_file):
        """
        Extract field of view from FITS WCS solution.
        
        Args:
            fits_file: Path to FITS file with WCS
            
        Returns:
            FOV in arcminutes or None
        """
        try:
            with fits.open(fits_file) as hdul:
                header = hdul[0].header
                
                # Try to get pixel scale from CD matrix or CDELT
                if 'CD1_1' in header and 'CD2_2' in header:
                    pixel_scale = np.sqrt(header['CD1_1']**2 + header['CD1_2']**2) * 3600  # arcsec/pixel
                elif 'CDELT1' in header:
                    pixel_scale = abs(header['CDELT1']) * 3600  # arcsec/pixel
                else:
                    return None
                
                # Get image dimensions
                naxis1 = header.get('NAXIS1', 0)
                naxis2 = header.get('NAXIS2', 0)
                
                if naxis1 > 0 and naxis2 > 0:
                    # FOV is diagonal of the image
                    diagonal_pixels = np.sqrt(naxis1**2 + naxis2**2)
                    fov_arcsec = diagonal_pixels * pixel_scale
                    fov_arcmin = fov_arcsec / 60.0
                    return fov_arcmin
                    
        except Exception as e:
            print(f"    Error extracting FOV: {e}")
        
        return None
        
    def query_vsp_photometry(self, target_ra, target_dec, fov_arcmin, target_name=None):
        """
        Query AAVSO VSP for comparison stars.
        
        Args:
            target_ra: Target RA in degrees
            target_dec: Target Dec in degrees  
            fov_arcmin: Field of view in arcminutes
            target_name: Optional target name
            
        Returns:
            List of comparison stars or None if query fails
        """
        try:
            # VSP API endpoint - NOTE: Currently not working per user
            url = "https://www.aavso.org/apps/vsp/api/chart/"
            
            params = {
                'format': 'json',
                'fov': fov_arcmin,
                'maglimit': 17.0,  # Magnitude limit
                'ra': target_ra,
                'dec': target_dec
            }
            
            if target_name:
                params['star'] = target_name
            
            print(f"    Querying VSP for comparison stars...")
            print(f"    Note: VSP API may not be working - attempting anyway")
            
            response = requests.get(url, params=params, timeout=30)
            
            if response.status_code == 200:
                data = response.json()
                
                if 'photometry' in data:
                    comp_stars = []
                    for star in data['photometry']:
                        if star.get('comp', False):  # Only comparison stars
                            comp_stars.append({
                                'auid': star.get('auid'),
                                'ra': star.get('ra'),
                                'dec': star.get('dec'),
                                'mag_v': star.get('v_mag'),
                                'mag_b': star.get('b_mag'),
                                'label': star.get('label')
                            })
                    
                    print(f"    Found {len(comp_stars)} comparison stars")
                    return comp_stars
            
            print(f"    VSP query failed (status {response.status_code})")
            return None
            
        except Exception as e:
            print(f"    Error querying VSP: {e}")
            return None
    
    def process_science_targets(self):
        """Process all unprocessed science targets (Phase 3)."""
        cursor = self.conn.cursor()
        
        # Get all targets that need processing (proc_bias and proc_flat are 0)
        cursor.execute('''
            SELECT date, target_name, filter, num_images 
            FROM target 
            WHERE proc_bias = 0 AND proc_flat = 0
            ORDER BY date, target_name, filter
        ''')
        
        targets = cursor.fetchall()
        
        if not targets:
            print("\nNo targets need processing")
            return
            
        print(f"\n{'='*60}")
        print(f"PHASE 3: Processing {len(targets)} target(s)")
        print(f"{'='*60}")
        
        for date, target_name, filter_name, num_images in targets:
            print(f"\nProcessing {date}/{target_name} (filter: {filter_name}, {num_images} images)")
            
            # Find calibration files
            bias_date, bias_path = self.find_closest_calibration(date, 'Bias')
            flat_date, flat_path = self.find_closest_calibration(date, 'Flat', filter_name)
            
            if bias_path:
                print(f"  Using bias from {bias_date}: {bias_path}")
            else:
                print(f"  WARNING: No bias calibration found")
                bias_date = 'nocal'
                
            if flat_path:
                print(f"  Using flat from {flat_date}: {flat_path}")
            else:
                print(f"  WARNING: No flat calibration found")
                flat_date = 'nocal'
            
            # Get science files
            light_dir = self.data_dir / date / 'Light' / target_name
            if not light_dir.exists():
                print(f"  ERROR: Light directory not found: {light_dir}")
                continue
                
            # Create Reduced_images directory
            reduced_dir = light_dir / 'Reduced_images'
            reduced_dir.mkdir(exist_ok=True)
            
            science_files = sorted(light_dir.glob('*.fits'))
            # Filter by the correct filter
            science_files = [f for f in science_files if self.get_filter_from_fits(f) == filter_name]
            
            processed_count = 0
            for science_file in science_files:
                # Skip if already in Reduced_images
                if 'Reduced_images' in str(science_file):
                    continue
                    
                output_file = reduced_dir / science_file.name
                
                # Calibrate the frame
                print(f"  Calibrating: {science_file.name}")
                if self.calibrate_frame(science_file, bias_path, flat_path, output_file):
                    # Plate solve
                    print(f"  Plate solving: {science_file.name}")
                    solution = self.plate_solve_astap(output_file)
                    
                    # Record WCS solution
                    cursor.execute('''
                        INSERT OR REPLACE INTO wcs 
                        (date, target_name, filter, filename, success, ra_deg, dec_deg, 
                         solve_time, num_stars, faintest_mag, fov_deg, error_message, timestamp)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    ''', (
                        date, target_name, filter_name, science_file.name,
                        1 if solution['success'] else 0,
                        solution.get('ra_deg'),
                        solution.get('dec_deg'),
                        solution.get('solve_time'),
                        solution.get('num_stars'),
                        solution.get('faintest_mag'),
                        solution.get('fov_deg'),
                        solution.get('error_message'),
                        solution['timestamp']
                    ))
                    
                    if solution['success']:
                        print(f"    ✓ Solved: RA={solution.get('ra_deg', 0):.4f}°, Dec={solution.get('dec_deg', 0):.4f}°")
                    else:
                        print(f"    ✗ Failed: {solution.get('error_message', 'Unknown error')}")
                        
                    processed_count += 1
                else:
                    print(f"    ✗ Calibration failed")
            
            # Update target table
            cursor.execute('''
                UPDATE target 
                SET proc_bias = ?, proc_flat = ?
                WHERE date = ? AND target_name = ? AND filter = ?
            ''', (bias_date, flat_date, date, target_name, filter_name))
            
            self.conn.commit()
            print(f"  Processed {processed_count}/{len(science_files)} frames")
    
    def process_phase4_photometry(self):
        """Process Phase 4: VSX queries and comparison star retrieval."""
        cursor = self.conn.cursor()
        
        # Get all targets that have been plate solved (have WCS data)
        cursor.execute('''
            SELECT DISTINCT date, target_name, filter
            FROM wcs
            WHERE success = 1
            ORDER BY date, target_name, filter
        ''')
        
        targets = cursor.fetchall()
        
        if not targets:
            print("\nNo plate-solved targets found for Phase 4 processing")
            return
            
        print(f"\n{'='*60}")
        print(f"PHASE 4: Processing {len(targets)} target(s) for photometry")
        print(f"{'='*60}")
        
        for date, target_name, filter_name in targets:
            print(f"\nPhase 4: {date}/{target_name} (filter: {filter_name})")
            
            # Check if already processed
            cursor.execute('''
                SELECT apphot FROM target
                WHERE date = ? AND target_name = ? AND filter = ?
            ''', (date, target_name, filter_name))
            
            result = cursor.fetchone()
            if result and result[0] not in [None, '', 0, 'nocal']:
                print(f"  Already processed (apphot={result[0]})")
                continue
            
            # Query VSX for target information
            print(f"  Querying VSX for {target_name}...")
            vsx_info = self.query_vsx(target_name)
            
            if not vsx_info:
                print(f"  ✗ Target not found in VSX")
                cursor.execute('''
                    UPDATE target SET apphot = 'nocomps'
                    WHERE date = ? AND target_name = ? AND filter = ?
                ''', (date, target_name, filter_name))
                self.conn.commit()
                continue
            
            print(f"  ✓ VSX: {vsx_info['name']} (AUID: {vsx_info['auid']})")
            print(f"    Type: {vsx_info['var_type']}, Period: {vsx_info['period']} d")
            print(f"    Mag: {vsx_info['mag_max']} - {vsx_info['mag_min']}")
            print(f"    Coordinates: RA={vsx_info['ra_deg']:.4f}°, Dec={vsx_info['dec_deg']:.4f}°")
            
            # Get a representative FITS file for this target/filter
            light_dir = self.data_dir / date / 'Light' / target_name / 'Reduced_images'
            if not light_dir.exists():
                print(f"  ✗ Reduced images directory not found")
                continue
            
            fits_files = [f for f in light_dir.glob('*.fits') if self.get_filter_from_fits(f) == filter_name]
            if not fits_files:
                print(f"  ✗ No FITS files found")
                continue
            
            # Use first file for calculations
            sample_fits = fits_files[0]
            
            # Calculate time and airmass for the sample
            time_info = self.calculate_time_and_airmass(
                sample_fits, 
                vsx_info['ra_deg'], 
                vsx_info['dec_deg']
            )
            
            if time_info:
                print(f"  Time: JD={time_info['jd']:.5f}, HJD={time_info['hjd']:.5f}")
                print(f"  Airmass: {time_info['airmass']:.3f}, Alt={time_info['altitude']:.1f}°")
            
            # Get FOV from WCS
            fov = self.get_fov_from_wcs(sample_fits)
            if fov:
                print(f"  FOV: {fov:.1f} arcmin")
            else:
                print(f"  Warning: Could not determine FOV, using default 60 arcmin")
                fov = 60.0
            
            # Query VSP for comparison stars
            comp_stars = self.query_vsp_photometry(
                vsx_info['ra_deg'],
                vsx_info['dec_deg'],
                fov,
                target_name
            )
            
            if comp_stars and len(comp_stars) > 0:
                print(f"  ✓ VSP returned {len(comp_stars)} comparison stars")
                apphot_status = 'ready'
            else:
                print(f"  ✗ No comparison stars available (VSP may be down)")
                apphot_status = 'nocomps'
            
            # Update target table
            cursor.execute('''
                UPDATE target SET apphot = ?
                WHERE date = ? AND target_name = ? AND filter = ?
            ''', (apphot_status, date, target_name, filter_name))
            
            self.conn.commit()
            print(f"  Updated apphot status: {apphot_status}")
        
    def run(self, phase3_only=False, phase4_only=False):
        """
        Main pipeline execution.
        
        Args:
            phase3_only: If True, skip phases 1 & 2 and only process science targets
            phase4_only: If True, skip phases 1-3 and only process photometry setup
        """
        print(f"Starting Wombat Pipeline")
        print(f"Data directory: {self.data_dir}")
        
        if not self.data_dir.exists():
            print(f"Error: Data directory {self.data_dir} does not exist")
            return
            
        self.setup_database()
        
        if phase4_only:
            # Skip everything and just do Phase 4
            self.process_phase4_photometry()
        elif phase3_only:
            # Skip phases 1 & 2, do phases 3 & 4
            self.process_science_targets()
            self.process_phase4_photometry()
        else:
            # Full pipeline: all phases
            # Get all night directories (subdirectories in data/)
            night_dirs = [d for d in self.data_dir.iterdir() if d.is_dir()]
            night_dirs.sort()
            
            print(f"\nFound {len(night_dirs)} night directories")
            
            for night_dir in night_dirs:
                self.process_night(night_dir)
                
            print("\n" + "="*60)
            print("Phase 1 & 2: Directory scan and calibration complete!")
            print("="*60)
            
            # Phase 3: Process science targets
            self.process_science_targets()
            
            # Phase 4: Photometry setup
            self.process_phase4_photometry()
        
        print("\n" + "="*60)
        print("Pipeline processing complete!")
        print("="*60)
        
        # Print summary
        self.print_summary()
        
    def print_summary(self):
        """Print summary of database contents."""
        cursor = self.conn.cursor()
        
        print("\nTarget Summary:")
        cursor.execute('SELECT date, target_name, filter, num_images FROM target ORDER BY date, target_name, filter')
        targets = cursor.fetchall()
        for date, target, filt, num_imgs in targets:
            print(f"  {date} - {target} ({filt}): {num_imgs} images")
            
        print("\nCalibration Summary:")
        cursor.execute('SELECT date, cal_type, filter, num_images FROM cal ORDER BY date, cal_type, filter')
        cals = cursor.fetchall()
        for date, cal_type, filt, num in cals:
            if filt:
                print(f"  {date} - {cal_type} ({filt}): {num} images")
            else:
                print(f"  {date} - {cal_type}: {num} images")
                
    def close(self):
        """Close database connection."""
        if self.conn:
            self.conn.close()


if __name__ == '__main__':
    pipeline = WombatPipeline()
    try:
        pipeline.run()
    finally:
        pipeline.close()
