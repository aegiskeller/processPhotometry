#!/usr/bin/env python3
"""
Web interface for Wombat Pipeline - Validation & Publishing
"""

from flask import Flask, render_template, jsonify, request, send_file
from pathlib import Path
import sqlite3
import json
from datetime import datetime
from datetime import datetime

app = Flask(__name__)
app.config['TEMPLATES_AUTO_RELOAD'] = True
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0

# Database path
DB_PATH = Path(__file__).parent / "wombatpipeline.db"


def get_db_connection():
    """Create database connection."""
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn


@app.route('/')
def index():
    """Database views page - Targets, WCS, Calibrations."""
    return render_template('database.html')


@app.route('/validate')
def validate():
    """Validate & Publish tab."""
    return render_template('validate.html')


@app.route('/api/targets/validate')
def get_validate_targets():
    """Get list of targets with apphot='validate' status."""
    conn = get_db_connection()
    cursor = conn.cursor()
    
    cursor.execute('''
        SELECT date, target_name, filter, num_images, submitted
        FROM target
        WHERE apphot = 'validate'
        ORDER BY date DESC, target_name
    ''')
    
    targets = []
    for row in cursor.fetchall():
        targets.append({
            'date': row['date'],
            'target_name': row['target_name'],
            'filter': row['filter'],
            'num_images': row['num_images'],
            'submitted': row['submitted']
        })
    
    conn.close()
    return jsonify(targets)


@app.route('/api/lightcurve/<date>/<target_name>/<filter_name>')
def get_lightcurve_data(date, target_name, filter_name):
    """Get lightcurve data for plotting."""
    
    # Try to find the photometry table in the standard directory structure
    # Format: D:\{date}\Light\{target_name}\Reduced_images\{target_name}_{filter}_photometry_transformed.tbl
    
    # Try multiple possible locations
    data_roots = [
        Path('D:/'),
        Path('C:/Users/Admin/Pictures/SharpCap Captures'),
        Path.cwd()
    ]
    
    phot_file = None
    for root in data_roots:
        # Try with _transformed suffix
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
        return jsonify({'error': f'Photometry file not found for {target_name} on {date}'}), 404
    
    # Get AUID from database
    conn = get_db_connection()
    cursor = conn.cursor()
    cursor.execute('''
        SELECT chart_id FROM target
        WHERE date = ? AND target_name = ? AND filter = ?
    ''', (date, target_name, filter_name))
    result = cursor.fetchone()
    auid = result['chart_id'] if result and result['chart_id'] else 'UNKNOWN'
    conn.close()
    
    # Read ZP log file to get rejection status
    zp_log_file = phot_file.parent / f'{target_name}_{filter_name}_photometry_zp.log'
    rejected_indices = set()
    
    if zp_log_file.exists():
        try:
            with open(zp_log_file, 'r') as f:
                for line in f:
                    if line.startswith('#') or line.startswith('Frame') or line.startswith('---'):
                        continue
                    if 'REJECTED' in line:
                        # Extract frame number from filename
                        # Format: 2025-09-30_19-49-47_V_-10.00_30.00s_0043.fits
                        parts = line.split()
                        if parts:
                            filename = parts[0]
                            # Extract frame number from filename (last part before .fits)
                            if '_' in filename:
                                try:
                                    frame_num = int(filename.split('_')[-1].replace('.fits', ''))
                                    rejected_indices.add(frame_num)
                                except (ValueError, IndexError):
                                    pass
        except Exception as e:
            print(f"Warning: Could not read ZP log file: {e}")
    
    try:
        # Read photometry table
        with open(phot_file, 'r') as f:
            lines = f.readlines()
        
        if len(lines) < 2:
            return jsonify({'error': 'Empty photometry file'}), 400
        
        # Parse header
        header = lines[0].strip().split('\t')
        col_idx = {col: i for i, col in enumerate(header)}
        
        # Check for transformed or instrumental magnitudes
        has_transformed = 'V_std_T1' in col_idx and 'V_std_T2' in col_idx
        
        if has_transformed:
            t1_mag_col = 'V_std_T1'
            t2_mag_col = 'V_std_T2'
            mag_type = 'Standard'
        else:
            t1_mag_col = 'Source_AMag_T1'
            t2_mag_col = 'Source_AMag_T2'
            mag_type = 'Instrumental'
        
        # Required columns
        required = ['HJD_UTC', t1_mag_col, 'Source_AMag_Err_T1', t2_mag_col, 'Source_AMag_Err_T2']
        missing = [col for col in required if col not in col_idx]
        if missing:
            return jsonify({'error': f'Missing columns: {missing}'}), 400
        
        # Extract data
        t1_data = []
        t2_data = []
        
        for i, line in enumerate(lines[1:]):
            if not line.strip():
                continue
            
            values = line.strip().split('\t')
            
            try:
                hjd = float(values[col_idx['HJD_UTC']])
                
                # T1 data
                t1_mag = float(values[col_idx[t1_mag_col]])
                t1_err = float(values[col_idx['Source_AMag_Err_T1']])
                
                t1_data.append({
                    'index': i,
                    'hjd': hjd,
                    'mag': t1_mag,
                    'err': t1_err,
                    'rejected': i in rejected_indices
                })
                
                # T2 data
                t2_mag = float(values[col_idx[t2_mag_col]])
                t2_err = float(values[col_idx['Source_AMag_Err_T2']])
                
                t2_data.append({
                    'index': i,
                    'hjd': hjd,
                    'mag': t2_mag,
                    'err': t2_err,
                    'rejected': i in rejected_indices
                })
                
            except (ValueError, IndexError) as e:
                continue
        
        return jsonify({
            'target_name': target_name,
            'auid': auid,
            'date': date,
            'filter': filter_name,
            'mag_type': mag_type,
            'T1': t1_data,
            'T2': t2_data
        })
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/target/update_status', methods=['POST'])
def update_target_status():
    """Update target submission status."""
    data = request.json
    
    date = data.get('date')
    target_name = data.get('target_name')
    filter_name = data.get('filter')
    status = data.get('status')  # 'usrrej', 'approved', 'submitted'
    
    if not all([date, target_name, filter_name, status]):
        return jsonify({'error': 'Missing parameters'}), 400
    
    conn = get_db_connection()
    cursor = conn.cursor()
    
    try:
        # Update both submitted status and apphot when approving
        if status == 'approved':
            cursor.execute('''
                UPDATE target SET submitted = ?, apphot = ?
                WHERE date = ? AND target_name = ? AND filter = ?
            ''', (status, 'approved', date, target_name, filter_name))
        else:
            cursor.execute('''
                UPDATE target SET submitted = ?
                WHERE date = ? AND target_name = ? AND filter = ?
            ''', (status, date, target_name, filter_name))
        
        conn.commit()
        
        if cursor.rowcount == 0:
            return jsonify({'error': 'Target not found'}), 404
        
        return jsonify({'success': True})
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500
    finally:
        conn.close()


@app.route('/api/target/save_rejected_points', methods=['POST'])
def save_rejected_points():
    """Save rejected points for a target."""
    data = request.json
    
    date = data.get('date')
    target_name = data.get('target_name')
    filter_name = data.get('filter')
    rejected_indices = data.get('rejected_indices', [])
    
    # Save to a JSON file for this target
    rejection_file = Path(f'rejections_{target_name}_{date}_{filter_name}.json')
    
    with open(rejection_file, 'w') as f:
        json.dump({
            'date': date,
            'target_name': target_name,
            'filter': filter_name,
            'rejected_indices': rejected_indices,
            'timestamp': datetime.now().isoformat()
        }, f, indent=2)
    
    return jsonify({'success': True, 'file': str(rejection_file)})


@app.route('/api/aavso/rebuild', methods=['POST'])
def rebuild_aavso_report():
    """Rebuild AAVSO report excluding rejected points."""
    data = request.json
    
    date = data.get('date')
    target_name = data.get('target_name')
    filter_name = data.get('filter')
    rejected_indices = set(data.get('rejected_indices', []))
    
    # Find the photometry table using same logic as lightcurve endpoint
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
        return jsonify({'error': f'Photometry file not found for {target_name} on {date}'}), 404
    
    # Get variable AUID and chart ID from database
    conn = get_db_connection()
    cursor = conn.cursor()
    cursor.execute('''
        SELECT variable_auid, chart_id FROM target
        WHERE date = ? AND target_name = ? AND filter = ?
    ''', (date, target_name, filter_name))
    result = cursor.fetchone()
    variable_auid = result['variable_auid'] if result and result['variable_auid'] else 'UNKNOWN'
    chart_id = result['chart_id'] if result and result['chart_id'] else 'UNKNOWN'
    conn.close()
    
    # Read ZP log file to get sigma-clipped rejection status
    zp_log_file = phot_file.parent / f'{target_name}_{filter_name}_photometry_zp.log'
    sigma_rejected_indices = set()
    
    if zp_log_file.exists():
        try:
            with open(zp_log_file, 'r') as f:
                for line in f:
                    if line.startswith('#') or line.startswith('Frame') or line.startswith('---'):
                        continue
                    if 'REJECTED' in line:
                        # Extract frame number from filename
                        parts = line.split()
                        if parts:
                            filename = parts[0]
                            if '_' in filename:
                                try:
                                    frame_num = int(filename.split('_')[-1].replace('.fits', ''))
                                    sigma_rejected_indices.add(frame_num)
                                except (ValueError, IndexError):
                                    pass
        except Exception as e:
            print(f"Warning: Could not read ZP log file: {e}")
    
    try:
        # Read photometry table
        with open(phot_file, 'r') as f:
            lines = f.readlines()
        
        header = lines[0].strip().split('\t')
        col_idx = {col: i for i, col in enumerate(header)}
        
        # Check for transformed magnitudes
        has_transformed = 'V_std_T1' in col_idx
        
        if not has_transformed:
            return jsonify({'error': 'Transformed magnitudes not available'}), 400
        
        # Build approved report with only non-rejected points
        output_file = Path(f'AAVSO_Report_{target_name}_{date}_approved.txt')
        
        # Get check star info - need V_std_T2 column for check star magnitude
        check_auid = 'na'  # Default if not available
        
        with open(output_file, 'w') as out:
            # Write AAVSO header - use proper line breaks, not escaped \n
            out.write("#TYPE=Extended\n")
            out.write("#OBSCODE=KSCA\n")
            out.write("#SOFTWARE=Wombat Pipeline v1.0\n")
            out.write("#DELIM=,\n")
            out.write("#DATE=HJD\n")
            out.write("#OBSTYPE=CCD\n")
            out.write("#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES\n")
            
            # Write data for non-rejected points
            points_written = 0
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
                    
                    # Get airmass - it might be stored in wrong column, default to 1.0 if unreasonable
                    try:
                        airmass_raw = float(values[col_idx['AIRMASS']]) if 'AIRMASS' in col_idx else 1.0
                        # Sanity check - airmass should be between 1.0 and 3.0 typically
                        airmass = airmass_raw if 1.0 <= airmass_raw <= 10.0 else 1.0
                    except (ValueError, IndexError):
                        airmass = 1.0
                    
                    # Format: NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES
                    # Use proper line breaks, not escaped \n
                    # NAME should be variable AUID, CHART should be chart ID
                    out.write(f"{variable_auid},{hjd:.6f},{v_std_t1:.3f},{v_err_t1:.3f},V,YES,STD,ENSEMBLE,na,{check_auid},{v_std_t2:.3f},{airmass:.2f},1,{chart_id},\n")
                    points_written += 1
                    
                except (ValueError, IndexError, KeyError):
                    continue
        
        total_rejected = len(sigma_rejected_indices) + len(rejected_indices)
        return jsonify({
            'success': True, 
            'file': str(output_file), 
            'points_included': points_written,
            'sigma_rejected': len(sigma_rejected_indices),
            'user_rejected': len(rejected_indices),
            'total_rejected': total_rejected
        })
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/aavso/auth')
def aavso_auth():
    """Initiate AAVSO OAuth flow."""
    state = secrets.token_urlsafe(32)
    session['oauth_state'] = state
    
    # Store return URL (which target to submit after auth)
    session['submit_after_auth'] = request.args.get('target')
    session['submit_date'] = request.args.get('date')
    session['submit_filter'] = request.args.get('filter')
    
    auth_url = f"{AAVSO_AUTH_URL}?response_type=code&client_id={AAVSO_CLIENT_ID}&redirect_uri={AAVSO_REDIRECT_URI}&state={state}"
    return redirect(auth_url)

@app.route('/aavso/callback')
def aavso_callback():
    """Handle AAVSO OAuth callback."""
    # Verify state to prevent CSRF
    if request.args.get('state') != session.get('oauth_state'):
        return jsonify({'error': 'Invalid state parameter'}), 400
    
    code = request.args.get('code')
    if not code:
        return jsonify({'error': 'No authorization code received'}), 400
    
    # Exchange code for access token
    token_data = {
        'grant_type': 'authorization_code',
        'code': code,
        'redirect_uri': AAVSO_REDIRECT_URI,
        'client_id': AAVSO_CLIENT_ID
    }
    
    try:
        response = requests.post(AAVSO_TOKEN_URL, data=token_data)
        response.raise_for_status()
        token_info = response.json()
        
        # Store access token in session
        session['aavso_access_token'] = token_info.get('access_token')
        session['aavso_token_type'] = token_info.get('token_type', 'Bearer')
        
        # If there was a pending submission, redirect to submit it
        if session.get('submit_after_auth'):
            return redirect(url_for('submit_to_aavso', 
                                   target=session.get('submit_after_auth'),
                                   date=session.get('submit_date'),
                                   filter=session.get('submit_filter')))
        
        return redirect('/validate?auth=success')
    
    except Exception as e:
        return jsonify({'error': f'Token exchange failed: {str(e)}'}), 500

@app.route('/api/aavso/submit', methods=['POST'])
def submit_to_aavso():
    """Generate and return AAVSO report for submission."""
    data = request.json
    
    date = data.get('date')
    target_name = data.get('target_name')
    filter_name = data.get('filter')
    
    # Look for approved report first, then fall back to original
    report_files = list(Path('.').glob(f'AAVSO_AID_Report_*_{date}_approved.txt'))
    
    if not report_files:
        report_files = list(Path('.').glob(f'AAVSO_AID_Report_*_{date}.txt'))
    
    if not report_files:
        return jsonify({'error': 'AAVSO report not found'}), 404
    
    # Return the most recent report
    latest_report = sorted(report_files, key=lambda p: p.stat().st_mtime, reverse=True)[0]
    
    return send_file(latest_report, as_attachment=True)


@app.route('/api/targets/all')
def get_all_targets():
    """Get all targets from database"""
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    
    cursor.execute("""
        SELECT date, target_name, filter, num_images, proc_bias, proc_flat, apphot, submitted
        FROM target
        ORDER BY date DESC, target_name, filter
    """)
    
    targets = [dict(row) for row in cursor.fetchall()]
    conn.close()
    
    return jsonify(targets)


@app.route('/api/wcs/all')
def get_all_wcs():
    """Get all WCS solutions from database"""
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    
    cursor.execute("""
        SELECT date, target_name, filter, filename, success, ra_deg, dec_deg,
               solve_time, num_stars, faintest_mag, fov_deg, error_message, timestamp
        FROM wcs
        ORDER BY date DESC, target_name, filename
    """)
    
    wcs_solutions = [dict(row) for row in cursor.fetchall()]
    conn.close()
    
    return jsonify(wcs_solutions)


@app.route('/api/wcs/target/<date>/<target_name>/<filter>')
def get_wcs_for_target(date, target_name, filter):
    """Get WCS solutions for a specific target"""
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    
    cursor.execute("""
        SELECT date, target_name, filter, filename, success, ra_deg, dec_deg,
               solve_time, num_stars, faintest_mag, fov_deg, error_message, timestamp
        FROM wcs
        WHERE date = ? AND target_name = ? AND filter = ?
        ORDER BY filename
    """, (date, target_name, filter))
    
    wcs_solutions = [dict(row) for row in cursor.fetchall()]
    conn.close()
    
    return jsonify(wcs_solutions)


@app.route('/api/calibrations/all')
def get_all_calibrations():
    """Get all calibration frames from database"""
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    
    cursor.execute("""
        SELECT date, cal_type, filter, num_images
        FROM cal
        ORDER BY date DESC, cal_type, filter
    """)
    
    calibrations = [dict(row) for row in cursor.fetchall()]
    conn.close()
    
    return jsonify(calibrations)


@app.route('/api/stats')
def get_stats():
    """Get pipeline statistics"""
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    
    # Total targets
    cursor.execute("SELECT COUNT(*) FROM target")
    total_targets = cursor.fetchone()[0]
    
    # Total nights (distinct dates)
    cursor.execute("SELECT COUNT(DISTINCT date) FROM target")
    total_nights = cursor.fetchone()[0]
    
    conn.close()
    
    return jsonify({
        'totalTargets': total_targets,
        'totalNights': total_nights
    })


if __name__ == '__main__':
    print("Starting Wombat Pipeline Web Interface")
    print("=" * 60)
    print(f"Database: {DB_PATH}")
    print("Navigate to: http://localhost:5000")
    print("=" * 60)
    app.run(host='127.0.0.1', port=5000, debug=False)
