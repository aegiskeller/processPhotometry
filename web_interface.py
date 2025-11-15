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
    # Find the transformed photometry table
    # Expected format: {AUID}_{DATE}_phot_transformed.tbl
    
    photometry_files = list(Path('.').glob(f'*_{date}_phot_transformed.tbl'))
    
    if not photometry_files:
        # Try without _transformed
        photometry_files = list(Path('.').glob(f'*_{date}_phot.tbl'))
    
    if not photometry_files:
        return jsonify({'error': 'Photometry file not found'}), 404
    
    phot_file = photometry_files[0]
    
    # Extract AUID from filename (format: {AUID}_{DATE}_phot*.tbl)
    auid = phot_file.stem.split('_')[0]
    
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
                    'rejected': False
                })
                
                # T2 data
                t2_mag = float(values[col_idx[t2_mag_col]])
                t2_err = float(values[col_idx['Source_AMag_Err_T2']])
                
                t2_data.append({
                    'index': i,
                    'hjd': hjd,
                    'mag': t2_mag,
                    'err': t2_err,
                    'rejected': False
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
    
    # Find the photometry table
    photometry_files = list(Path('.').glob(f'*_{date}_phot_transformed.tbl'))
    if not photometry_files:
        photometry_files = list(Path('.').glob(f'*_{date}_phot.tbl'))
    
    if not photometry_files:
        return jsonify({'error': 'Photometry file not found'}), 404
    
    phot_file = photometry_files[0]
    auid = phot_file.stem.split('_')[0]
    
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
        output_file = Path(f'AAVSO_AID_Report_{auid}_{date}_approved.txt')
        
        with open(output_file, 'w') as out:
            # Write AAVSO header
            out.write("#TYPE=Extended\\n")
            out.write("#OBSCODE=YOUR_CODE\\n")
            out.write("#SOFTWARE=Wombat Pipeline v1.0\\n")
            out.write("#DELIM=,\\n")
            out.write("#DATE=JD\\n")
            out.write("#OBSTYPE=CCD\\n")
            out.write("#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES\\n")
            
            # Write data for non-rejected points
            for i, line in enumerate(lines[1:]):
                if not line.strip():
                    continue
                
                # Check if this index is rejected
                if f'T1_{i}' in rejected_indices:
                    continue
                
                values = line.strip().split('\t')
                
                try:
                    hjd = float(values[col_idx['HJD_UTC']])
                    v_std = float(values[col_idx['V_std_T1']])
                    v_err = float(values[col_idx['Source_AMag_Err_T1']])
                    airmass = float(values[col_idx['AIRMASS']]) if 'AIRMASS' in col_idx else 1.0
                    
                    # Format: NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES
                    out.write(f"{auid},{hjd:.6f},{v_std:.3f},{v_err:.3f},V,YES,STD,ENSEMBLE,na,na,na,{airmass:.2f},na,na,\\n")
                    
                except (ValueError, IndexError, KeyError):
                    continue
        
        return jsonify({'success': True, 'file': str(output_file), 'points_included': i + 1 - len(rejected_indices)})
        
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
    print("Navigate to: http://localhost:8080")
    print("=" * 60)
    app.run(host='0.0.0.0', port=8080, debug=False)
