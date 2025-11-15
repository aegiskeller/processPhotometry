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
from astropy.wcs import WCS
import astropy.units as u
import numpy as np
import subprocess
import re
from datetime import datetime
import requests
from urllib.parse import urlencode
from photutils.detection import DAOStarFinder
from photutils.background import MMMBackground, MADStdBackgroundRMS
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
from scipy.optimize import curve_fit


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
                        # VSX returns RA/Dec as decimal degrees
                        try:
                            ra_deg = float(ra_str)
                            dec_deg = float(dec_str)
                        except (ValueError, TypeError):
                            print(f"    Warning: Could not parse coordinates for {obj.get('Name')}")
                            return 'novsx'
                        
                        # Convert to HMS/DMS for display
                        coord = SkyCoord(ra=ra_deg*u.degree, dec=dec_deg*u.degree)
                        ra_hms = coord.ra.to_string(unit=u.hour, sep=':')
                        dec_dms = coord.dec.to_string(unit=u.degree, sep=':')
                        
                        # Parse magnitude values (format: "12.31 V" or "12.8" or "(0.3)")
                        mag_max_str = obj.get('MaxMag', '')
                        mag_min_str = obj.get('MinMag', '')
                        
                        # Extract numeric value from magnitude string
                        mag_max = None
                        mag_min = None
                        
                        try:
                            mag_min = float(mag_min_str.split()[0]) if mag_min_str else None
                        except (ValueError, IndexError):
                            mag_min = None
                        
                        # Handle MaxMag: could be absolute value or relative "(0.3)"
                        if mag_max_str:
                            if mag_max_str.strip().startswith('(') and mag_max_str.strip().endswith(')'):
                                # Relative amplitude: MaxMag = MinMag + amplitude
                                try:
                                    amplitude = float(mag_max_str.strip('()').split()[0])
                                    if mag_min is not None:
                                        mag_max = mag_min + amplitude
                                except (ValueError, IndexError):
                                    mag_max = None
                            else:
                                # Absolute magnitude value
                                try:
                                    mag_max = float(mag_max_str.split()[0])
                                except (ValueError, IndexError):
                                    mag_max = None
                        
                        return {
                            'name': obj.get('Name'),
                            'auid': obj.get('AUID'),
                            'ra_deg': ra_deg,
                            'dec_deg': dec_deg,
                            'ra_hms': ra_hms,
                            'dec_dms': dec_dms,
                            'var_type': obj.get('Type'),
                            'period': obj.get('Period'),
                            'mag_max': mag_max,
                            'mag_min': mag_min,
                            'mag_max_str': mag_max_str,
                            'mag_min_str': mag_min_str,
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
            # VSP API endpoint
            url = "https://www.aavso.org/apps/vsp/api/chart/"
            
            params = {
                'format': 'json',
                'fov': fov_arcmin,
                'maglimit': 17.0  # Magnitude limit
            }
            
            # Prefer star name if provided (more reliable than coordinates)
            if target_name:
                params['star'] = target_name
            else:
                params['ra'] = target_ra
                params['dec'] = target_dec
            
            print(f"    Querying VSP for comparison stars...")
            
            response = requests.get(url, params=params, timeout=30)
            
            if response.status_code == 200:
                data = response.json()
                
                if 'photometry' in data:
                    comp_stars = []
                    
                    # All stars in the photometry array are comparison stars
                    for star in data['photometry']:
                        # Parse RA/Dec from sexagesimal to decimal degrees
                        ra_str = star.get('ra', '')
                        dec_str = star.get('dec', '')
                        
                        try:
                            # Parse RA (HH:MM:SS.SS) to degrees
                            ra_parts = ra_str.split(':')
                            ra_decimal = (float(ra_parts[0]) + float(ra_parts[1]) / 60 + float(ra_parts[2]) / 3600) * 15
                            
                            # Parse Dec (DD:MM:SS.S)
                            dec_parts = dec_str.split(':')
                            dec_sign = -1 if dec_parts[0].startswith('-') else 1
                            dec_decimal = dec_sign * (abs(float(dec_parts[0])) + float(dec_parts[1]) / 60 + float(dec_parts[2]) / 3600)
                        except (ValueError, IndexError):
                            print(f"      Warning: Could not parse coordinates for star {star.get('auid')}")
                            continue
                        
                        comp_star = {
                            'auid': star.get('auid'),
                            'label': star.get('label'),
                            'ra': ra_decimal,
                            'dec': dec_decimal,
                            'chart_id': data.get('chartid')
                        }
                        
                        # Extract magnitudes from bands array
                        # VSP structure: 'bands': [{'band': 'V', 'mag': 11.908, 'error': 0.029}, ...]
                        for band_data in star.get('bands', []):
                            band = band_data.get('band')
                            mag = band_data.get('mag')
                            err = band_data.get('error')
                            
                            if band and mag is not None:
                                # Map VSP band names to our convention
                                # V, B, Rc, Ic → mag_v, mag_b, mag_r, mag_i
                                band_map = {'V': 'v', 'B': 'b', 'Rc': 'r', 'Ic': 'i'}
                                our_band = band_map.get(band, band.lower())
                                
                                comp_star[f'mag_{our_band}'] = mag
                                if err is not None:
                                    comp_star[f'mag_{our_band}_error'] = err
                        
                        comp_stars.append(comp_star)
                    
                    print(f"    ✓ Found {len(comp_stars)} comparison stars from VSP chart {data.get('chartid')}")
                    return comp_stars
            
            print(f"    ✗ VSP query failed (status {response.status_code})")
            return 'novsp'  # Service unavailable
            
        except Exception as e:
            print(f"    ✗ Error querying VSP: {e}")
            return 'novsp'  # Network error or service down
    
    def filter_comparison_stars(self, comp_stars, fits_file, vsx_info, filter_name='V'):
        """
        Filter comparison stars by image bounds and magnitude range (Phase 5).
        
        Args:
            comp_stars: List of comparison stars from VSP
            fits_file: Path to FITS file (for WCS and image dimensions)
            vsx_info: VSX target information (with mag_min, mag_max)
            filter_name: Filter to use for magnitude comparison (V, B, R, I)
            
        Returns:
            Filtered list of comparison stars, or empty list if none suitable
        """
        if not comp_stars or comp_stars == 'novsp':
            return []
        
        try:
            # Open FITS file to get WCS and image dimensions
            with fits.open(fits_file) as hdul:
                header = hdul[0].header
                
                # Get image dimensions
                naxis1 = header.get('NAXIS1')  # Width
                naxis2 = header.get('NAXIS2')  # Height
                
                if not naxis1 or not naxis2:
                    print("    Warning: Cannot determine image dimensions")
                    return []
                
                # Try to create WCS object
                try:
                    from astropy.wcs import WCS
                    wcs = WCS(header)
                except Exception as e:
                    print(f"    Warning: Cannot create WCS from header: {e}")
                    return []
            
            # Magnitude range: 1.5 mag brighter than min (brightest), 1.5 mag fainter than max (faintest)
            # Variable mag range is max (bright) to min (faint), so:
            if vsx_info.get('mag_max') is not None and vsx_info.get('mag_min') is not None:
                mag_bright_limit = vsx_info['mag_max'] - 1.5  # Don't use stars brighter than this
                mag_faint_limit = vsx_info['mag_min'] + 1.5   # Don't use stars fainter than this
                print(f"    Filtering comp stars: mag range {mag_bright_limit:.1f} - {mag_faint_limit:.1f}")
            else:
                # No magnitude limits available - use wide range
                mag_bright_limit = 8.0
                mag_faint_limit = 18.0
                print(f"    Warning: No VSX magnitude limits, using default range {mag_bright_limit:.1f} - {mag_faint_limit:.1f}")
            
            filtered_stars = []
            mag_key = f'mag_{filter_name.lower()}'
            
            for star in comp_stars:
                # Check if star has magnitude in this filter
                star_mag = star.get(mag_key)
                if star_mag is None:
                    continue
                
                # Filter by magnitude range
                if star_mag < mag_bright_limit or star_mag > mag_faint_limit:
                    continue
                
                # Convert RA/Dec to pixel coordinates
                try:
                    ra = star['ra']
                    dec = star['dec']
                    
                    # WCS world_to_pixel expects (ra, dec) in degrees
                    pixel_coords = wcs.world_to_pixel_values(ra, dec)
                    x, y = pixel_coords[0], pixel_coords[1]
                    
                    # Check if within image bounds (with small margin)
                    margin = 50  # pixels from edge
                    if (x < margin or x > naxis1 - margin or 
                        y < margin or y > naxis2 - margin):
                        continue
                    
                    # Star passed all filters
                    star['pixel_x'] = x
                    star['pixel_y'] = y
                    filtered_stars.append(star)
                    
                except Exception as e:
                    print(f"      Warning: Could not convert coords for star {star.get('auid')}: {e}")
                    continue
            
            print(f"    ✓ {len(filtered_stars)}/{len(comp_stars)} comp stars passed filters")
            return filtered_stars
            
        except Exception as e:
            print(f"    ✗ Error filtering comparison stars: {e}")
            return []
    
    def select_check_star(self, comp_stars, vsx_info, filter_name='V'):
        """
        Select check star (T2) from comparison stars based on proximity to variable.
        The check star should be closest in both position and magnitude to the variable.
        
        Args:
            comp_stars: List of filtered comparison stars (with pixel_x, pixel_y)
            vsx_info: VSX target information (with ra_deg, dec_deg, mag_max, mag_min)
            filter_name: Filter to use for magnitude comparison
            
        Returns:
            Dict with 'check_star' (T2) and 'comp_stars' (C3..N), or None if no suitable stars
        """
        if not comp_stars or len(comp_stars) == 0:
            return None
        
        try:
            # Get variable star properties
            var_ra = vsx_info.get('ra_deg', vsx_info.get('ra', 0))
            var_dec = vsx_info.get('dec_deg', vsx_info.get('dec', 0))
            var_mag = (vsx_info.get('mag_max', 12.0) + vsx_info.get('mag_min', 12.0)) / 2.0  # Mean magnitude
            
            mag_key = f'mag_{filter_name.lower()}'
            
            # Calculate combined distance metric for each comp star
            # Distance metric = spatial_distance (arcsec) + magnitude_difference * scaling_factor
            # Scale magnitude difference to be comparable to arcsec (1 mag ~ 60 arcsec)
            best_star = None
            best_metric = float('inf')
            
            for star in comp_stars:
                # Skip if no magnitude in this filter
                star_mag = star.get(mag_key)
                if star_mag is None:
                    continue
                
                # Calculate spatial separation in arcsec
                from astropy.coordinates import SkyCoord
                import astropy.units as u
                
                var_coord = SkyCoord(ra=var_ra*u.degree, dec=var_dec*u.degree)
                star_coord = SkyCoord(ra=star['ra']*u.degree, dec=star['dec']*u.degree)
                spatial_sep = var_coord.separation(star_coord).arcsecond
                
                # Calculate magnitude difference
                mag_diff = abs(star_mag - var_mag)
                
                # Combined metric: weight spatial and magnitude equally
                # Use scaling: 1 mag difference ~ 60 arcsec separation
                metric = spatial_sep + (mag_diff * 60.0)
                
                if metric < best_metric:
                    best_metric = metric
                    best_star = star
                    best_spatial_sep = spatial_sep
            
            if best_star is None:
                return None
            
            # Prepare result
            check_star = best_star.copy()
            check_star['label'] = 'T2'
            
            # Label remaining stars as C3, C4, ...
            comp_list = []
            comp_index = 3
            for star in comp_stars:
                if star['auid'] != check_star['auid']:
                    labeled_star = star.copy()
                    labeled_star['label'] = f'C{comp_index}'
                    comp_list.append(labeled_star)
                    comp_index += 1
            
            print(f"    ✓ Selected check star T2: {check_star['auid']}")
            print(f"      Position: RA={check_star['ra']:.5f}, Dec={check_star['dec']:.5f}")
            print(f"      Magnitude: {check_star.get(mag_key, 'N/A')}")
            print(f"      Separation from variable: {best_spatial_sep:.1f} arcsec")
            print(f"    ✓ {len(comp_list)} additional comp stars labeled C3-C{comp_index-1}")
            
            return {
                'check_star': check_star,
                'comp_stars': comp_list
            }
            
        except Exception as e:
            print(f"    ✗ Error selecting check star: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def perform_aperture_photometry(self, fits_files, stars_dict, vsx_info, aperture_params, filter_name='V', output_file=None):
        """
        Perform aperture photometry on all stars across all images.
        
        Args:
            fits_files: List of FITS file paths
            stars_dict: Dict with 'check_star' (T2), 'comp_stars' (C3-CN)
            vsx_info: VSX target information
            aperture_params: Dict with aperture_radius, annulus_inner, annulus_outer, fwhm
            filter_name: Filter used (for output naming)
            output_file: Optional output file path (default: auto-generated)
            
        Returns:
            Path to output .tbl file or None on error
        """
        try:
            # Prepare star list: T1 (variable), T2 (check), C3-CN (comps)
            all_stars = []
            
            # T1: Variable star
            var_star = {
                'label': 'T1',
                'ra': vsx_info.get('ra_deg'),
                'dec': vsx_info.get('dec_deg'),
                'auid': vsx_info.get('auid'),
                'name': vsx_info.get('name')
            }
            all_stars.append(var_star)
            
            # T2: Check star
            check_star = stars_dict['check_star'].copy()
            all_stars.append(check_star)
            
            # C3-CN: Comparison stars
            for comp in stars_dict['comp_stars']:
                all_stars.append(comp)
            
            print(f"\n    Performing aperture photometry on {len(all_stars)} stars across {len(fits_files)} images...")
            print(f"    Apertures: r={aperture_params['aperture_radius']:.0f} px, annulus: {aperture_params['annulus_inner']:.0f}-{aperture_params['annulus_outer']:.0f} px")
            
            # Observatory location (need to get from FITS or config)
            # Using a default location - should be configurable
            observatory = EarthLocation.of_site('Kitt Peak')
            
            # Photometry results for all images
            phot_table_rows = []
            
            for img_num, fits_file in enumerate(sorted(fits_files), 1):
                print(f"    Processing {img_num}/{len(fits_files)}: {fits_file.name}")
                
                try:
                    with fits.open(fits_file) as hdul:
                        data = hdul[0].data.astype(float)
                        header = hdul[0].header
                        wcs = WCS(header)
                        
                        # Get observation time and calculate JD, HJD, BJD
                        obs_time_str = header.get('DATE-OBS') or header.get('DATE-AVG')
                        if obs_time_str:
                            obs_time = Time(obs_time_str, format='isot', scale='utc')
                        else:
                            obs_time = Time(header.get('MJD-OBS', 0), format='mjd', scale='utc')
                        
                        jd_utc = obs_time.jd
                        
                        # Calculate HJD and BJD for target position
                        target_coord = SkyCoord(ra=vsx_info['ra_deg']*u.deg, dec=vsx_info['dec_deg']*u.deg)
                        ltt_bary = obs_time.light_travel_time(target_coord, 'barycentric', location=observatory)
                        ltt_helio = obs_time.light_travel_time(target_coord, 'heliocentric', location=observatory)
                        
                        bjd_tdb = (obs_time.tdb + ltt_bary).jd
                        hjd_utc = (obs_time.utc + ltt_helio).jd
                        
                        # Calculate airmass
                        altaz = target_coord.transform_to(AltAz(obstime=obs_time, location=observatory))
                        airmass = altaz.secz.value
                        altitude = altaz.alt.deg
                        
                        # Get exposure time and other metadata
                        exptime = header.get('EXPTIME', header.get('EXPOSURE', 0))
                        ccd_temp = header.get('CCD-TEMP', header.get('SET-TEMP', -999))
                        
                        # Convert star RA/Dec to pixel positions and perform photometry
                        row_data = {
                            'Label': fits_file.name,
                            'slice': img_num,
                            'Saturated': 0,  # Could check for saturation
                            'J.D.-2400000': jd_utc - 2400000,
                            'JD_UTC': jd_utc,
                            'JD_SOBS': jd_utc,  # Start of observation
                            'HJD_UTC': hjd_utc,
                            'BJD_TDB': bjd_tdb,
                            'AIRMASS': airmass,
                            'ALT_OBJ': altitude,
                            'CCD-TEMP': ccd_temp,
                            'EXPTIME': exptime,
                            'EXPOSURE': exptime,
                            'RAOBJ2K': vsx_info['ra_deg'] / 15.0,  # Convert to hours
                            'DECOBJ2K': vsx_info['dec_deg'],
                            'FWHM_Mean': aperture_params['fwhm'],
                            'Source_Radius': aperture_params['aperture_radius'],
                            'Sky_Rad(min)': aperture_params['annulus_inner'],
                            'Sky_Rad(max)': aperture_params['annulus_outer']
                        }
                        
                        # Perform photometry on each star
                        for star in all_stars:
                            label = star['label']
                            
                            # Convert RA/Dec to pixel coordinates
                            pixel_coords = wcs.world_to_pixel_values(star['ra'], star['dec'])
                            x, y = pixel_coords[0], pixel_coords[1]
                            
                            # Create apertures
                            position = [(x, y)]
                            aperture = CircularAperture(position, r=aperture_params['aperture_radius'])
                            annulus = CircularAnnulus(position, 
                                                     r_in=aperture_params['annulus_inner'],
                                                     r_out=aperture_params['annulus_outer'])
                            
                            # Perform aperture photometry
                            phot = aperture_photometry(data, aperture)
                            annulus_phot = aperture_photometry(data, annulus)
                            
                            # Calculate sky background per pixel
                            sky_per_pixel = annulus_phot['aperture_sum'][0] / annulus.area
                            
                            # Sky-subtracted source flux
                            source_minus_sky = phot['aperture_sum'][0] - (sky_per_pixel * aperture.area)
                            
                            # Calculate error (Poisson + readnoise)
                            # Simplified error model - should include gain, readnoise from header
                            gain = header.get('GAIN', 1.0)
                            readnoise = header.get('RDNOISE', 5.0)
                            
                            # Error in source flux
                            n_pix_source = aperture.area
                            n_pix_sky = annulus.area
                            
                            source_error = np.sqrt(
                                phot['aperture_sum'][0] / gain +  # Poisson noise in source
                                n_pix_source * readnoise**2 +      # Readnoise in source aperture
                                (n_pix_source / n_pix_sky) * annulus_phot['aperture_sum'][0] / gain +  # Poisson in sky
                                (n_pix_source / n_pix_sky) * readnoise**2  # Readnoise in sky measurement
                            ) * gain
                            
                            # Instrumental magnitude
                            if source_minus_sky > 0:
                                source_amag = -2.5 * np.log10(source_minus_sky)
                                source_amag_err = 1.0857 * source_error / source_minus_sky  # 1.0857 = 2.5/ln(10)
                                snr = source_minus_sky / source_error
                            else:
                                source_amag = 99.999
                                source_amag_err = 9.999
                                snr = 0.0
                            
                            # Get peak value in source aperture
                            cutout_size = int(aperture_params['aperture_radius']) + 2
                            x_int, y_int = int(x), int(y)
                            if (y_int - cutout_size >= 0 and y_int + cutout_size < data.shape[0] and
                                x_int - cutout_size >= 0 and x_int + cutout_size < data.shape[1]):
                                cutout = data[y_int-cutout_size:y_int+cutout_size+1, 
                                             x_int-cutout_size:x_int+cutout_size+1]
                                peak_value = np.max(cutout)
                                mean_value = np.mean(cutout)
                            else:
                                peak_value = 0
                                mean_value = 0
                            
                            # Store measurements for this star
                            row_data[f'X(IJ)_{label}'] = x
                            row_data[f'Y(IJ)_{label}'] = y
                            row_data[f'X(FITS)_{label}'] = x + 1  # FITS is 1-indexed
                            row_data[f'Y(FITS)_{label}'] = y + 1
                            row_data[f'RA_{label}'] = star['ra'] / 15.0  # Convert to hours
                            row_data[f'DEC_{label}'] = star['dec']
                            row_data[f'Source-Sky_{label}'] = source_minus_sky
                            row_data[f'Source_Error_{label}'] = source_error
                            row_data[f'Source_AMag_{label}'] = source_amag
                            row_data[f'Source_AMag_Err_{label}'] = source_amag_err
                            row_data[f'Source_SNR_{label}'] = snr
                            row_data[f'Peak_{label}'] = peak_value
                            row_data[f'Mean_{label}'] = mean_value
                            row_data[f'Sky/Pixel_{label}'] = sky_per_pixel
                            row_data[f'N_Src_Pixels_{label}'] = int(n_pix_source)
                            row_data[f'N_Sky_Pixels_{label}'] = int(n_pix_sky)
                            row_data[f'FWHM_{label}'] = aperture_params['fwhm']
                            row_data[f'Width_{label}'] = aperture_params['fwhm']
                            row_data[f'X-Width_{label}'] = aperture_params['fwhm']
                            row_data[f'Y-Width_{label}'] = aperture_params['fwhm']
                            row_data[f'Angle_{label}'] = 0.0
                            row_data[f'Roundness_{label}'] = 1.0
                        
                        # Calculate relative fluxes (normalized to sum of comparison stars)
                        comp_labels = [s['label'] for s in all_stars if s['label'].startswith('C')]
                        total_comp_flux = sum(row_data[f'Source-Sky_C{i}'] for i in range(3, 3+len(comp_labels)))
                        total_comp_error = np.sqrt(sum(row_data[f'Source_Error_C{i}']**2 for i in range(3, 3+len(comp_labels))))
                        
                        row_data['tot_C_cnts'] = total_comp_flux
                        row_data['tot_C_err'] = total_comp_error
                        
                        for star in all_stars:
                            label = star['label']
                            source_flux = row_data[f'Source-Sky_{label}']
                            source_err = row_data[f'Source_Error_{label}']
                            
                            if total_comp_flux > 0:
                                rel_flux = source_flux / total_comp_flux
                                rel_flux_err = rel_flux * np.sqrt(
                                    (source_err / source_flux)**2 + (total_comp_error / total_comp_flux)**2
                                )
                                rel_flux_snr = rel_flux / rel_flux_err if rel_flux_err > 0 else 0
                            else:
                                rel_flux = 0
                                rel_flux_err = 0
                                rel_flux_snr = 0
                            
                            row_data[f'rel_flux_{label}'] = rel_flux
                            row_data[f'rel_flux_err_{label}'] = rel_flux_err
                            row_data[f'rel_flux_SNR_{label}'] = rel_flux_snr
                        
                        phot_table_rows.append(row_data)
                        
                except Exception as e:
                    print(f"      Warning: Failed to process {fits_file.name}: {e}")
                    import traceback
                    traceback.print_exc()
                    continue
            
            if not phot_table_rows:
                print(f"    ✗ No images successfully processed")
                return None
            
            print(f"    ✓ Aperture photometry complete for {len(phot_table_rows)} images")
            
            # Generate output table file
            if output_file is None:
                # Auto-generate filename: AUID_DATE_phot.tbl
                date_str = datetime.now().strftime('%Y%m%d')
                auid = vsx_info.get('auid', 'unknown').replace('-', '')
                output_file = Path(f"{auid}_{date_str}_phot.tbl")
            
            # Write photometry table in TSV format
            self._write_photometry_table(output_file, phot_table_rows, all_stars)
            
            print(f"    ✓ Photometry table written to: {output_file}")
            return output_file
            
        except Exception as e:
            print(f"    ✗ Error in aperture photometry: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def create_aavso_report(self, photometry_table_file, obscode, chart_id=None, notes="", date_type='HJD', obstype='CCD'):
        """
        Create AAVSO Variable Star Report from photometry table.
        
        Args:
            photometry_table_file: Path to .tbl photometry file
            obscode: AAVSO observer code (e.g., 'KSCA')
            chart_id: VSP chart identifier (from query_vsp_photometry)
            notes: Observing conditions notes
            date_type: 'JD' or 'HJD' (default: 'HJD')
            obstype: 'CCD' or 'DSLR' (default: 'CCD')
            
        Returns:
            Path to AAVSO report file or None on error
        """
        try:
            from datetime import datetime
            
            print(f"\n    Creating AAVSO Variable Star Report...")
            print(f"    Observer Code: {obscode}")
            print(f"    Date Type: {date_type}")
            
            # Read photometry table
            with open(photometry_table_file, 'r') as f:
                lines = f.readlines()
            
            if len(lines) < 2:
                print(f"    ✗ Empty or invalid photometry table")
                return None
            
            # Parse header
            header = lines[0].strip().split('\t')
            
            # Create column index map
            col_idx = {col: i for i, col in enumerate(header)}
            
            # Extract metadata from table
            filter_used = 'V'  # Default, could extract from filename or table
            target_name = None
            
            # Try to extract target name from filename (format: AUID_DATE_phot.tbl)
            table_path = Path(photometry_table_file)
            name_parts = table_path.stem.split('_')
            if len(name_parts) >= 2:
                # First part might be AUID or target name
                target_name = name_parts[0]
            
            if not target_name:
                print(f"    ✗ Cannot determine target name from file")
                return None
            
            # AAVSO Extended Format specification
            TYPE = "Extended"
            SOFTWARE = f"AstroImageJ Python Pipeline"
            DELIM = ","
            
            # Determine date column based on date_type
            if date_type == 'HJD':
                date_col = 'HJD_UTC'
            elif date_type == 'BJD':
                date_col = 'BJD_TDB'
            else:  # JD
                date_col = 'JD_UTC'
            
            # Verify required columns exist
            required_cols = [date_col, 'AIRMASS', 'EXPTIME', 'Source_AMag_T1', 'Source_AMag_Err_T1']
            missing_cols = [col for col in required_cols if col not in col_idx]
            if missing_cols:
                print(f"    ✗ Missing required columns: {missing_cols}")
                return None
            
            # Determine check star (T2) and comparison stars
            # Check star is T2
            check_star_amag_col = 'Source_AMag_T2'
            check_star_sky_col = 'Source-Sky_T2'
            
            # Find comparison stars (C3, C4, C5, ...)
            comp_star_cols = [col for col in header if col.startswith('Source-Sky_C')]
            n_comp_stars = len(comp_star_cols)
            
            print(f"    Found {n_comp_stars} comparison stars")
            
            # Check star information (from VSP query - would be stored)
            # For now, use placeholder values
            check_star_label = "T2"
            check_star_cat_mag = "NA"
            check_star_cat_mag_err = "NA"
            
            # Determine if using ensemble or single comparison star
            if n_comp_stars == 1:
                # Single comparison star
                comp_label = comp_star_cols[0].replace('Source-Sky_', '')
                use_ensemble = False
            else:
                # Ensemble of comparison stars
                use_ensemble = True
                comp_label = "ENSEMBLE"
            
            # Generate output filename
            timestamp = datetime.now().strftime('%Y%m%d-%Hh%Mm%Ss')
            output_file = Path(f"AAVSO_AID_Report_{obscode}_{target_name}_{date_type}_{timestamp}.txt")
            
            # Write AAVSO report
            with open(output_file, 'w') as f:
                # Write header keywords
                f.write(f"#TYPE={TYPE}\n")
                f.write(f"#OBSCODE={obscode}\n")
                f.write(f"#SOFTWARE={SOFTWARE}\n")
                f.write(f"#DELIM={DELIM}\n")
                f.write(f"#DATE={date_type}\n")
                f.write(f"#OBSTYPE={obstype}\n")
                
                # Write data header
                data_header = "#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES"
                f.write(data_header + "\n")
                
                # Write data rows
                n_rows = 0
                for line in lines[1:]:  # Skip header
                    if not line.strip():
                        continue
                    
                    values = line.strip().split('\t')
                    
                    # DATE: JD/HJD/BJD with appropriate precision
                    date_value = f"{float(values[col_idx[date_col]]):.6f}"
                    
                    # MAG: Target apparent magnitude (3 decimal places)
                    mag = f"{float(values[col_idx['Source_AMag_T1']]):.3f}"
                    
                    # MERR: Target magnitude error (3 decimal places)
                    merr = f"{float(values[col_idx['Source_AMag_Err_T1']]):.3f}"
                    
                    # FILT: Filter used
                    filt = filter_used
                    
                    # TRANS: Transformation applied (NO for untransformed)
                    trans = "NO"
                    
                    # MTYPE: Magnitude type (STD = standard, DIF = differential)
                    mtype = "STD"
                    
                    # CNAME: Comparison star name or ENSEMBLE
                    if use_ensemble:
                        cname = "ENSEMBLE"
                        cmag = "na"
                    else:
                        cname = comp_label
                        # Instrumental magnitude of comp star
                        comp_sky = float(values[col_idx[comp_star_cols[0]]])
                        exptime = float(values[col_idx['EXPTIME']])
                        cmag_value = -2.5 * np.log10(comp_sky / exptime)
                        cmag = f"{cmag_value:.3f}"
                    
                    # KNAME: Check star name
                    kname = check_star_label
                    
                    # KMAG: Check star magnitude
                    if use_ensemble:
                        # Use apparent magnitude
                        kmag = f"{float(values[col_idx[check_star_amag_col]]):.3f}"
                    else:
                        # Use instrumental magnitude
                        check_sky = float(values[col_idx[check_star_sky_col]])
                        exptime = float(values[col_idx['EXPTIME']])
                        kmag_value = -2.5 * np.log10(check_sky / exptime)
                        kmag = f"{kmag_value:.3f}"
                    
                    # AMASS: Airmass
                    amass = f"{float(values[col_idx['AIRMASS']]):.6f}"
                    
                    # GROUP: Grouping identifier (use 1 for all)
                    group = "1"
                    
                    # CHART: VSP chart identifier
                    chart = chart_id if chart_id else "NA"
                    
                    # NOTES: Extended information
                    check_amag = f"{float(values[col_idx[check_star_amag_col]]):.3f}"
                    check_sky = float(values[col_idx[check_star_sky_col]])
                    exptime = float(values[col_idx['EXPTIME']])
                    check_ins_mag = -2.5 * np.log10(check_sky / exptime)
                    
                    target_sky = float(values[col_idx['Source-Sky_T1']])
                    target_ins_mag = -2.5 * np.log10(target_sky / exptime)
                    
                    notes_parts = [
                        notes,
                        f"KMAG={check_amag}",
                        f"KMAGINS={check_ins_mag:.3f}",
                        f"KREFMAG={check_star_cat_mag}",
                        f"KREFERR={check_star_cat_mag_err}",
                        f"VMAGINS={target_ins_mag:.3f}"
                    ]
                    
                    # If single comp star, add CMAGINS
                    if not use_ensemble:
                        notes_parts.append(f"CMAGINS={cmag}")
                    
                    notes_str = "|".join([n for n in notes_parts if n])
                    
                    # Construct data row
                    data_row = f"{target_name},{date_value},{mag},{merr},{filt},{trans},{mtype},{cname},{cmag},{kname},{kmag},{amass},{group},{chart},{notes_str}"
                    f.write(data_row + "\n")
                    n_rows += 1
            
            print(f"    ✓ AAVSO report created: {output_file}")
            print(f"    {n_rows} observations written")
            
            return output_file
            
        except Exception as e:
            print(f"    ✗ Error creating AAVSO report: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def transform_to_standard_magnitudes(self, photometry_table_file, vsp_matches, color_basis='vi'):
        """
        Transform instrumental magnitudes to standard V magnitudes using comparison stars.
        
        Uses the transformation equation:
        V_standard = V_instrumental + ZP + k * color
        
        Where:
        - ZP is the zero point (offset)
        - k is the color coefficient (slope)
        - color is the color index (B-V, V-I, or V-R)
        
        Args:
            photometry_table_file: Path to .tbl photometry file
            vsp_matches: Dict mapping star labels (T2, C3, C4...) to VSP star data
                        Each entry should have 'bands' with magnitudes
            color_basis: Color index to use ('bv', 'vi', 'vr')
            
        Returns:
            Dict with transformation results or None on error
        """
        try:
            import statistics
            
            print(f"\n    Computing photometric transformation...")
            print(f"    Color basis: {color_basis.upper()}")
            
            # Read photometry table
            with open(photometry_table_file, 'r') as f:
                lines = f.readlines()
            
            if len(lines) < 2:
                print(f"    ✗ Empty or invalid photometry table")
                return None
            
            # Parse header and rows
            header = lines[0].strip().split('\t')
            col_idx = {col: i for i, col in enumerate(header)}
            
            rows = []
            for line in lines[1:]:
                if line.strip():
                    values = line.strip().split('\t')
                    row = {header[i]: values[i] for i in range(len(values))}
                    rows.append(row)
            
            # Collect transformation data points from matched comparison stars
            data_points = []
            
            for row in rows:
                for label, aavso_data in vsp_matches.items():
                    # Skip T1 (variable) - we'll apply transformation to it later
                    if label == 'T1':
                        continue
                    
                    # Get instrumental magnitude column
                    inst_mag_col = f"Source_AMag_{label}"
                    if inst_mag_col not in col_idx:
                        continue
                    
                    try:
                        inst_mag = float(row[inst_mag_col])
                    except (ValueError, TypeError, KeyError):
                        continue
                    
                    # Get standard magnitudes from VSP data
                    bands = aavso_data.get('bands', {})
                    std_v_mag = bands.get('V')
                    std_b_mag = bands.get('B')
                    std_i_mag = bands.get('Ic') or bands.get('I')
                    std_r_mag = bands.get('R')
                    
                    if std_v_mag is None:
                        continue
                    
                    # Calculate color index based on color_basis
                    color_index = None
                    if color_basis == 'bv' and std_b_mag is not None:
                        color_index = std_b_mag - std_v_mag
                    elif color_basis == 'vi' and std_i_mag is not None:
                        color_index = std_v_mag - std_i_mag
                    elif color_basis == 'vr' and std_r_mag is not None:
                        color_index = std_v_mag - std_r_mag
                    
                    if color_index is not None:
                        # Delta mag = standard - instrumental
                        delta_mag = std_v_mag - inst_mag
                        data_points.append({
                            'label': label,
                            'inst_mag': inst_mag,
                            'std_mag': std_v_mag,
                            'color': color_index,
                            'delta': delta_mag
                        })
            
            if len(data_points) < 2:
                print(f"    ✗ Insufficient data points for transformation ({len(data_points)} found, need ≥2)")
                return None
            
            print(f"    Collected {len(data_points)} data points from comparison stars")
            
            # Fit linear transformation: delta_mag = ZP + k * color
            # Using least squares fit
            n = len(data_points)
            sum_color = sum(p['color'] for p in data_points)
            sum_delta = sum(p['delta'] for p in data_points)
            sum_color_sq = sum(p['color']**2 for p in data_points)
            sum_color_delta = sum(p['color'] * p['delta'] for p in data_points)
            
            # Solve for slope (k) and intercept (ZP)
            denom = n * sum_color_sq - sum_color**2
            if abs(denom) < 1e-10:
                # All colors are the same, just compute mean zero point
                zero_point = sum_delta / n
                color_coeff = 0.0
            else:
                color_coeff = (n * sum_color_delta - sum_color * sum_delta) / denom
                zero_point = (sum_delta - color_coeff * sum_color) / n
            
            # Calculate RMS residual
            residuals = []
            for p in data_points:
                predicted = zero_point + color_coeff * p['color']
                residuals.append(p['delta'] - predicted)
            
            rms = np.sqrt(sum(r**2 for r in residuals) / len(residuals)) if residuals else float('nan')
            
            # Count unique stars
            unique_stars = len(set(p['label'] for p in data_points))
            
            print(f"    ✓ Transformation computed:")
            print(f"      Zero point: {zero_point:.3f}")
            print(f"      Color coefficient (k): {color_coeff:.3f}")
            print(f"      RMS residual: {rms:.3f} mag")
            print(f"      Stars used: {unique_stars}")
            
            # Compute median color from comparison stars for applying to T1 and T2
            comp_colors = [p['color'] for p in data_points]
            median_color = statistics.median(comp_colors) if comp_colors else 0.0
            
            print(f"      Median comparison color: {median_color:.3f}")
            
            # Apply transformation to T1 (variable) and T2 (check star)
            transformed_mags = {'T1': [], 'T2': []}
            
            for row_idx, row in enumerate(rows):
                # Transform T1 using median comparison color
                try:
                    t1_inst = float(row['Source_AMag_T1'])
                    t1_std = t1_inst + zero_point + color_coeff * median_color
                    transformed_mags['T1'].append({
                        'row': row_idx,
                        'inst_mag': t1_inst,
                        'std_mag': t1_std,
                        'color_used': median_color
                    })
                except (KeyError, ValueError, TypeError):
                    pass
                
                # Transform T2 using median comparison color
                try:
                    t2_inst = float(row['Source_AMag_T2'])
                    t2_std = t2_inst + zero_point + color_coeff * median_color
                    transformed_mags['T2'].append({
                        'row': row_idx,
                        'inst_mag': t2_inst,
                        'std_mag': t2_std,
                        'color_used': median_color
                    })
                except (KeyError, ValueError, TypeError):
                    pass
            
            print(f"    ✓ Transformed {len(transformed_mags['T1'])} T1 measurements")
            print(f"    ✓ Transformed {len(transformed_mags['T2'])} T2 measurements")
            
            # Write transformed photometry table
            output_file = Path(str(photometry_table_file).replace('.tbl', '_transformed.tbl'))
            
            with open(output_file, 'w') as f:
                # Write header with additional columns
                new_header = header + ['V_std_T1', 'V_std_T2']
                f.write('\t'.join(new_header) + '\n')
                
                # Write rows with transformed magnitudes
                for row_idx, row in enumerate(rows):
                    # Get transformed mags for this row
                    t1_std = None
                    t2_std = None
                    
                    for t1_entry in transformed_mags['T1']:
                        if t1_entry['row'] == row_idx:
                            t1_std = t1_entry['std_mag']
                            break
                    
                    for t2_entry in transformed_mags['T2']:
                        if t2_entry['row'] == row_idx:
                            t2_std = t2_entry['std_mag']
                            break
                    
                    # Write original row values
                    row_values = [row.get(col, '') for col in header]
                    
                    # Add transformed magnitudes
                    row_values.append(f'{t1_std:.6f}' if t1_std is not None else '')
                    row_values.append(f'{t2_std:.6f}' if t2_std is not None else '')
                    
                    f.write('\t'.join(row_values) + '\n')
            
            print(f"    ✓ Transformed table written: {output_file}")
            
            return {
                'zero_point': zero_point,
                'color_coefficient': color_coeff,
                'rms': rms,
                'num_stars': unique_stars,
                'num_obs': len(data_points),
                'median_color': median_color,
                'transformed_mags': transformed_mags,
                'output_file': str(output_file),
                'success': True
            }
            
        except Exception as e:
            print(f"    ✗ Error in transformation: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def update_target_status(self, date, target_name, filter_name, apphot_status):
        """
        Update the apphot status for a target in the database.
        
        Args:
            date: Observation date (YYYYMMDD format)
            target_name: Name of the target
            filter_name: Filter used
            apphot_status: Status to set ('validate', 'nocomps', 'novsp', etc.)
        """
        try:
            cursor = self.conn.cursor()
            cursor.execute('''
                UPDATE target SET apphot = ?
                WHERE date = ? AND target_name = ? AND filter = ?
            ''', (apphot_status, date, target_name, filter_name))
            
            self.conn.commit()
            rows_updated = cursor.rowcount
            
            if rows_updated > 0:
                print(f"    ✓ Updated target status: apphot='{apphot_status}'")
                return True
            else:
                print(f"    ✗ No target found for {date}/{target_name}/{filter_name}")
                return False
                
        except Exception as e:
            print(f"    ✗ Error updating target status: {e}")
            return False
    
    def process_phase6_transformation(self, date, target_name, filter_name, photometry_table_file, vsp_matches, color_basis='vi'):
        """
        Complete Phase 6 workflow: transform magnitudes and update database status.
        
        Args:
            date: Observation date (YYYYMMDD format)
            target_name: Name of the target
            filter_name: Filter used
            photometry_table_file: Path to .tbl photometry file
            vsp_matches: Dict mapping star labels to VSP star data
            color_basis: Color index to use ('bv', 'vi', 'vr')
            
        Returns:
            True if successful and status updated to 'validate', False otherwise
        """
        try:
            print(f"\n{'='*60}")
            print(f"Phase 6: Photometric Transformation")
            print(f"Target: {date}/{target_name}/{filter_name}")
            print(f"{'='*60}")
            
            # Perform transformation
            result = self.transform_to_standard_magnitudes(
                photometry_table_file=photometry_table_file,
                vsp_matches=vsp_matches,
                color_basis=color_basis
            )
            
            if result and result.get('success'):
                # Check transformation quality
                rms = result.get('rms', float('inf'))
                num_stars = result.get('num_stars', 0)
                
                # Quality thresholds (adjust as needed)
                rms_threshold = 0.5  # magnitudes
                min_stars = 2
                
                if rms > rms_threshold:
                    print(f"    ⚠ Warning: RMS residual ({rms:.3f}) exceeds threshold ({rms_threshold})")
                    apphot_status = 'high_rms'
                elif num_stars < min_stars:
                    print(f"    ⚠ Warning: Too few stars ({num_stars}) for transformation")
                    apphot_status = 'few_stars'
                else:
                    print(f"    ✓ Transformation quality acceptable")
                    apphot_status = 'validate'
                
                # Update database
                success = self.update_target_status(date, target_name, filter_name, apphot_status)
                
                if success and apphot_status == 'validate':
                    print(f"\n{'='*60}")
                    print(f"✓ Phase 6 Complete: Status set to 'validate'")
                    print(f"{'='*60}")
                    return True
                else:
                    print(f"\n{'='*60}")
                    print(f"⚠ Phase 6 Complete: Status set to '{apphot_status}'")
                    print(f"{'='*60}")
                    return False
            else:
                # Transformation failed
                print(f"    ✗ Transformation failed")
                self.update_target_status(date, target_name, filter_name, 'trans_failed')
                return False
                
        except Exception as e:
            print(f"    ✗ Error in Phase 6 processing: {e}")
            import traceback
            traceback.print_exc()
            self.update_target_status(date, target_name, filter_name, 'trans_error')
            return False
    
    def _write_photometry_table(self, output_file, rows, stars):
        """Write photometry table in tab-separated format."""
        # Build column headers
        base_cols = ['Label', 'slice', 'Saturated', 'J.D.-2400000', 'JD_UTC', 'JD_SOBS', 
                     'HJD_UTC', 'BJD_TDB', 'AIRMASS', 'ALT_OBJ', 'CCD-TEMP', 'EXPTIME', 
                     'EXPOSURE', 'RAOBJ2K', 'DECOBJ2K', 'FWHM_Mean', 'Source_Radius', 
                     'Sky_Rad(min)', 'Sky_Rad(max)']
        
        # Relative flux columns for all stars
        rel_flux_cols = []
        for star in stars:
            label = star['label']
            rel_flux_cols.extend([f'rel_flux_{label}'])
        for star in stars:
            label = star['label']
            rel_flux_cols.extend([f'rel_flux_err_{label}'])
        for star in stars:
            label = star['label']
            rel_flux_cols.extend([f'rel_flux_SNR_{label}'])
        
        # Total comparison star flux
        comp_flux_cols = ['tot_C_cnts', 'tot_C_err']
        
        # Per-star detailed columns
        star_cols = []
        for star in stars:
            label = star['label']
            star_cols.extend([
                f'X(IJ)_{label}', f'Y(IJ)_{label}', f'X(FITS)_{label}', f'Y(FITS)_{label}',
                f'RA_{label}', f'DEC_{label}', f'Source-Sky_{label}', f'Source_Error_{label}',
                f'Source_AMag_{label}', f'Source_AMag_Err_{label}', f'Source_SNR_{label}',
                f'Peak_{label}', f'Mean_{label}', f'Sky/Pixel_{label}', f'FWHM_{label}',
                f'Width_{label}', f'X-Width_{label}', f'Y-Width_{label}', f'Angle_{label}',
                f'Roundness_{label}'
            ])
        
        # N_Src_Pixels and N_Sky_Pixels columns (only for T2 onwards per AIJ format)
        n_pix_cols = []
        for star in stars[1:]:  # Skip T1
            label = star['label']
            n_pix_cols.extend([f'N_Src_Pixels_{label}', f'N_Sky_Pixels_{label}'])
        
        all_cols = base_cols + rel_flux_cols + comp_flux_cols + star_cols + n_pix_cols
        
        # Write to file
        with open(output_file, 'w') as f:
            # Write header
            f.write('\t'.join(all_cols) + '\n')
            
            # Write data rows
            for row in rows:
                values = []
                for col in all_cols:
                    val = row.get(col, '')
                    if isinstance(val, float):
                        if 'JD' in col or 'AIRMASS' in col or 'ALT' in col:
                            values.append(f'{val:.10f}')
                        elif 'flux' in col or 'Source' in col or 'Sky' in col or 'Peak' in col:
                            values.append(f'{val:.10f}')
                        elif 'Mag' in col:
                            values.append(f'{val:.6f}')
                        elif 'SNR' in col:
                            values.append(f'{val:.10f}')
                        elif 'Width' in col or 'FWHM' in col:
                            values.append(f'{val:.10f}')
                        elif 'X(' in col or 'Y(' in col or 'RA_' in col or 'DEC_' in col:
                            values.append(f'{val:.10f}')
                        else:
                            values.append(f'{val:.10f}')
                    elif isinstance(val, int):
                        values.append(str(val))
                    else:
                        values.append(str(val))
                f.write('\t'.join(values) + '\n')
    
    def fit_psf_and_determine_aperture(self, fits_file, min_counts=5000, max_counts=30000, force_fwhm=None):
        """
        Find a suitable reference star and measure FWHM using radial profile (AIJ-style).
        
        Args:
            fits_file: Path to FITS file
            min_counts: Minimum peak counts for reference star
            max_counts: Maximum peak counts for reference star
            force_fwhm: If provided, use this FWHM value instead of measuring
            
        Returns:
            Dictionary with aperture parameters or None if measurement fails
        """
        try:
            if force_fwhm:
                print(f"    Using specified FWHM={force_fwhm:.2f} pixels")
                fwhm = force_fwhm
            else:
                print(f"    Determining aperture from seeing profile...")
            
            with fits.open(fits_file) as hdul:
                data = hdul[0].data.astype(float)
                
                # Estimate background
                bkg_estimator = MMMBackground()
                bkg = bkg_estimator(data)
                
                # Find stars in the image
                std_estimator = MADStdBackgroundRMS()
                std = std_estimator(data)
                
                daofind = DAOStarFinder(fwhm=4.0, threshold=5.0 * std)
                sources = daofind(data - bkg)
                
                if sources is None or len(sources) == 0:
                    print(f"    ✗ No stars detected in image")
                    return None
                
                # Find stars with peak counts in desired range
                suitable_stars = []
                for src in sources:
                    x, y = int(src['xcentroid']), int(src['ycentroid'])
                    
                    # Check if within image bounds
                    if x < 50 or x >= data.shape[1] - 50 or y < 50 or y >= data.shape[0] - 50:
                        continue
                    
                    # Get peak value in 5x5 region around centroid
                    cutout = data[y-2:y+3, x-2:x+3]
                    peak = np.max(cutout) - bkg
                    
                    if min_counts <= peak <= max_counts:
                        suitable_stars.append({
                            'x': x,
                            'y': y,
                            'peak': peak,
                            'flux': src['flux'],
                            'sharpness': src['sharpness']
                        })
                
                if not suitable_stars:
                    print(f"    ✗ No stars found with peak counts {min_counts}-{max_counts}")
                    return None
                
                # AIJ uses the brightest star with reasonable sharpness (0.3-0.7)
                # Filter for well-focused stars first
                well_focused = [s for s in suitable_stars if 0.3 <= s['sharpness'] <= 0.7]
                
                if well_focused:
                    # Use brightest well-focused star
                    ref_star = max(well_focused, key=lambda s: s['peak'])
                else:
                    # Fallback to brightest star regardless of sharpness
                    ref_star = max(suitable_stars, key=lambda s: s['peak'])
                
                print(f"    Reference star: ({ref_star['x']}, {ref_star['y']}), peak={ref_star['peak']:.0f} counts, sharpness={ref_star['sharpness']:.3f}")
                
                # Extract cutout around star
                size = 40
                x, y = ref_star['x'], ref_star['y']
                cutout = data[y-size:y+size+1, x-size:x+size+1].copy()
                cutout -= bkg
                
                # AIJ Seeing Profile method: Compute radial profile
                # Create radial distance array from center
                cy, cx = size, size  # Center of cutout
                y_idx, x_idx = np.ogrid[0:cutout.shape[0], 0:cutout.shape[1]]
                
                # Refine centroid: find actual peak within central region
                central = cutout[cy-2:cy+3, cx-2:cx+3]
                peak_offset = np.unravel_index(np.argmax(central), central.shape)
                cy_refined = cy + (peak_offset[0] - 2)
                cx_refined = cx + (peak_offset[1] - 2)
                
                # Recompute radii from refined center
                r = np.sqrt((x_idx - cx_refined)**2 + (y_idx - cy_refined)**2)
                
                # Compute radial bins (AIJ uses integer bins)
                max_radius = min(size, 30)
                radii_sum = np.zeros(max_radius)
                means_sum = np.zeros(max_radius)
                counts = np.zeros(max_radius, dtype=int)
                
                # Accumulate values in bins
                for j in range(cutout.shape[0]):
                    for i in range(cutout.shape[1]):
                        radius = r[j, i]
                        bin_idx = int(radius)
                        if bin_idx < max_radius:
                            radii_sum[bin_idx] += radius
                            means_sum[bin_idx] += cutout[j, i]
                            counts[bin_idx] += 1
                
                # Compute mean values per bin
                radii = np.zeros(max_radius)
                profile = np.zeros(max_radius)
                
                for bin_idx in range(max_radius):
                    if counts[bin_idx] > 0:
                        radii[bin_idx] = radii_sum[bin_idx] / counts[bin_idx]
                        profile[bin_idx] = means_sum[bin_idx] / counts[bin_idx]
                
                # Find mean peak (maximum of binned profile)
                mean_peak = np.max(profile)
                
                # Normalize profile by mean peak (AIJ method)
                profile_normalized = profile / mean_peak
                
                # Find FWHM: radius where normalized profile crosses 0.5
                try:
                    # Find where profile crosses half-max (0.5 after normalization)
                    above_half = profile_normalized > 0.5
                    if np.any(above_half):
                        # Find last index above half-max
                        last_above = np.where(above_half)[0][-1]
                        
                        if last_above < len(radii) - 1:
                            # Linear interpolation between points (AIJ method)
                            r1, r2 = radii[last_above], radii[last_above + 1]
                            v1, v2 = profile_normalized[last_above], profile_normalized[last_above + 1]
                            
                            # Interpolate to 0.5
                            if v1 != v2:
                                r_half = r1 + (0.5 - v1) * (r2 - r1) / (v2 - v1)
                            else:
                                r_half = r1
                            
                            fwhm = 2.0 * r_half
                        else:
                            # Use last point as estimate
                            fwhm = 2.0 * radii[last_above]
                    else:
                        # Peak not well-defined
                        print(f"    ✗ Could not measure FWHM from radial profile")
                        return 'badpsf'
                    
                    print(f"    Seeing profile: FWHM={fwhm:.2f} pixels")
                    
                    # Check if FWHM is reasonable (2-20 pixels typical range)
                    if fwhm < 2.0 or fwhm > 20.0:
                        print(f"    ✗ FWHM outside reasonable range (2-20 pixels)")
                        return 'badpsf'
                    
                    # AIJ aperture determination: find where flux drops below cutoff (default 1%)
                    # This is the PRIMARY method - FWHM fallback only if this fails
                    flux_cutoff = 0.01  # 1% of peak
                    found_r1 = False
                    
                    for bin_idx in range(1, len(profile_normalized)):
                        if profile_normalized[bin_idx] < flux_cutoff:
                            # Found where flux drops below 1%
                            r1 = radii[bin_idx]
                            r2 = r1 * 2.0
                            r3 = r2 * 2.0
                            found_r1 = True
                            break
                    
                    if not found_r1:
                        # Fallback to FWHM-based (AIJ lines 285-288)
                        r1 = fwhm * 1.7
                        r2 = fwhm * 3.4
                        r3 = fwhm * 6.8
                        print(f"    Using FWHM-based radii (fallback)")
                    else:
                        print(f"    Using flux-cutoff based radii (1% point)")
                    
                    # Round radii (AIJ does this)
                    aperture_radius = int(np.ceil(r1))
                    annulus_inner = int(np.ceil(r2))
                    annulus_outer = int(np.ceil(r3))
                    
                    print(f"    ✓ Aperture: r={aperture_radius} px, annulus: {annulus_inner}-{annulus_outer} px")
                    
                    return {
                        'fwhm': fwhm,
                        'aperture_radius': float(aperture_radius),
                        'annulus_inner': float(annulus_inner),
                        'annulus_outer': float(annulus_outer),
                        'ref_star_x': ref_star['x'],
                        'ref_star_y': ref_star['y']
                    }
                    
                except Exception as e:
                    print(f"    ✗ FWHM measurement failed: {e}")
                    return 'badpsf'
                    
        except Exception as e:
            print(f"    ✗ Error in seeing profile analysis: {e}")
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
            
            # Store directory name (may differ from AUID)
            dir_name = target_name
            
            # Update target_name with AUID if different from directory name
            auid = vsx_info['auid']
            if auid and auid != target_name:
                print(f"  Updating target_name from '{target_name}' to '{auid}'")
                cursor.execute('''
                    UPDATE target SET target_name = ?
                    WHERE date = ? AND target_name = ? AND filter = ?
                ''', (auid, date, target_name, filter_name))
                cursor.execute('''
                    UPDATE wcs SET target_name = ?
                    WHERE date = ? AND target_name = ? AND filter = ?
                ''', (auid, date, target_name, filter_name))
                self.conn.commit()
                target_name = auid  # Use updated name for database queries
            
            # Get a representative FITS file for this target/filter (using directory name)
            light_dir = self.data_dir / date / 'Light' / dir_name / 'Reduced_images'
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
            
            # Handle VSP service down
            if comp_stars == 'novsp':
                print(f"  ✗ VSP service unavailable")
                apphot_status = 'novsp'
            elif not comp_stars or len(comp_stars) == 0:
                print(f"  ✗ No comparison stars returned from VSP")
                apphot_status = 'nocomps'
            else:
                print(f"  ✓ VSP returned {len(comp_stars)} comparison stars")
                
                # Phase 5: Filter comparison stars by image bounds and magnitude
                filtered_stars = self.filter_comparison_stars(
                    comp_stars,
                    sample_fits,
                    vsx_info,
                    filter_name
                )
                
                if not filtered_stars or len(filtered_stars) == 0:
                    print(f"  ✗ No suitable comparison stars after filtering")
                    apphot_status = 'nocomps'
                else:
                    print(f"  ✓ {len(filtered_stars)} suitable comparison stars available")
                    apphot_status = 'ready'
            
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
