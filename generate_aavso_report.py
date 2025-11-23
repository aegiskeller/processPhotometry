#!/usr/bin/env python3
"""
Generate fresh AAVSO report with correct variable_auid and airmass values.
"""

import sys
sys.path.insert(0, '.')

from pathlib import Path
import sqlite3

# Import the rebuild function logic
date = '2025-09-30'
target_name = 'RW PsA'
filter_name = 'V'
rejected_indices = set()  # No user rejections

# Find the photometry table
data_roots = [
    Path('D:/'),
    Path('C:/Users/Admin/Pictures/SharpCap Captures'),
    Path.cwd()
]

phot_file = None
for root in data_roots:
    possible_paths = [
        root / date / 'Light' / target_name / 'Reduced_images' / f'{target_name}_{filter_name}_photometry_transformed.tbl',
        root / date / 'Light' / target_name / 'Reduced_images' / f'{target_name}_{filter_name}_photometry.tbl',
    ]
    
    for path in possible_paths:
        if path.exists():
            phot_file = path
            break
    
    if phot_file:
        break

if not phot_file:
    print(f'ERROR: Photometry file not found for {target_name} on {date}')
    sys.exit(1)

print(f"Using photometry file: {phot_file}")

# Get variable AUID and chart ID from database
conn = sqlite3.connect('wombatpipeline.db')
conn.row_factory = sqlite3.Row
cursor = conn.cursor()
cursor.execute('''
    SELECT variable_auid, chart_id FROM target
    WHERE date = ? AND target_name = ? AND filter = ?
''', (date, target_name, filter_name))
result = cursor.fetchone()
variable_auid = result['variable_auid'] if result and result['variable_auid'] else 'UNKNOWN'
chart_id = result['chart_id'] if result and result['chart_id'] else 'UNKNOWN'
conn.close()

print(f"Variable AUID: {variable_auid}")
print(f"Chart ID: {chart_id}")

# Read ZP log file to get sigma-clipped rejection status
zp_log_file = phot_file.parent / f'{target_name}_{filter_name}_photometry_zp.log'
sigma_rejected_indices = set()

if zp_log_file.exists():
    with open(zp_log_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('Frame') or line.startswith('---'):
                continue
            if 'REJECTED' in line:
                parts = line.split()
                if parts:
                    filename = parts[0]
                    if '_' in filename:
                        try:
                            frame_num = int(filename.split('_')[-1].replace('.fits', ''))
                            sigma_rejected_indices.add(frame_num)
                        except (ValueError, IndexError):
                            pass

print(f"Sigma-clipped frames: {len(sigma_rejected_indices)}")

# Read photometry table
with open(phot_file, 'r') as f:
    lines = f.readlines()

header = lines[0].strip().split('\t')
col_idx = {col: i for i, col in enumerate(header)}

# Check for transformed magnitudes
has_transformed = 'V_std_T1' in col_idx

if not has_transformed:
    print('ERROR: Transformed magnitudes not available')
    sys.exit(1)

# Build approved report with only non-rejected points
output_file = Path(f'AAVSO_Report_{target_name}_{date}_approved.txt')

check_auid = 'na'

points_written = 0
with open(output_file, 'w') as out:
    # Write AAVSO header
    out.write("#TYPE=Extended\n")
    out.write("#OBSCODE=KSCA\n")
    out.write("#SOFTWARE=Wombat Pipeline v1.0\n")
    out.write("#DELIM=,\n")
    out.write("#DATE=HJD\n")
    out.write("#OBSTYPE=CCD\n")
    out.write("#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES\n")
    
    # Write data for non-rejected points
    for i, line in enumerate(lines[1:]):
        if not line.strip():
            continue
        
        # Skip if rejected by sigma clipping
        if i in sigma_rejected_indices:
            continue
        
        # Skip if rejected by user in GUI
        if f'T1_{i}' in rejected_indices:
            continue
        
        values = line.strip().split('\t')
        
        try:
            hjd = float(values[col_idx['HJD_UTC']])
            v_std_t1 = float(values[col_idx['V_std_T1']])
            v_err_t1 = float(values[col_idx['Source_AMag_Err_T1']])
            v_std_t2 = float(values[col_idx['V_std_T2']]) if 'V_std_T2' in col_idx else 0.0
            
            # Get airmass
            try:
                airmass_raw = float(values[col_idx['AIRMASS']]) if 'AIRMASS' in col_idx else 1.0
                # Sanity check
                airmass = airmass_raw if 1.0 <= airmass_raw <= 10.0 else 1.0
            except (ValueError, IndexError):
                airmass = 1.0
            
            # Write line with variable_auid in NAME field and chart_id in CHART field
            out.write(f"{variable_auid},{hjd:.6f},{v_std_t1:.3f},{v_err_t1:.3f},V,YES,STD,ENSEMBLE,na,{check_auid},{v_std_t2:.3f},{airmass:.2f},1,{chart_id},\n")
            points_written += 1
            
        except (ValueError, IndexError, KeyError) as e:
            continue

total_rejected = len(sigma_rejected_indices) + len(rejected_indices)
print(f"\nReport generated: {output_file}")
print(f"Points included: {points_written}")
print(f"Sigma-clipped: {len(sigma_rejected_indices)}")
print(f"User rejected: {len(rejected_indices)}")
print(f"Total rejected: {total_rejected}")

# Display first 5 lines
print(f"\nFirst 5 data lines:")
with open(output_file, 'r') as f:
    lines = f.readlines()
    for i, line in enumerate(lines):
        if line.startswith('#'):
            continue
        if i < 12:  # Skip headers and show first 5 data lines
            print(line.rstrip())
