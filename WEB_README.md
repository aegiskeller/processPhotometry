# Wombat Pipeline - Validate & Publish Web Interface

Interactive web interface for validating photometry results and publishing to AAVSO.

## Features

### 🎯 Target Validation
- View all targets with `apphot='validate'` status
- Display target details: date, filter, number of images
- Track submission status (pending, approved, submitted, rejected)

### 📈 Interactive Lightcurves
- **T1 (Variable)**: Blue markers with error bars
- **T2 (Check star)**: Red markers with error bars
- **Inverted Y-axis**: Lower magnitudes at top (astronomical convention)
- Real-time plotting with Plotly.js

### 🔧 Interactive Tools
- **Zoom**: Scroll or box select
- **Pan**: Click and drag
- **Reset**: Double-click or use reset button
- **Point Selection**: Click individual points to toggle rejection
- **Box Select**: Drag to select multiple points for rejection
- **Lasso Select**: Available in mode bar

### ⚙️ Workflow Actions
1. **Select Target**: Click target from list to load lightcurve
2. **Review Data**: Examine T1 (variable) and T2 (check) lightcurves
3. **Reject Bad Points**: Click or box-select outliers
4. **Reject Target**: Mark entire dataset as rejected (`submitted='usrrej'`)
5. **Approve Target**: Mark for AAVSO submission (`submitted='approved'`)
6. **Push to AAVSO**: Download AAVSO report file (`submitted='submitted'`)

## Installation

### Prerequisites
```bash
pip install flask
```

Or install from requirements:
```bash
pip install -r requirements_web.txt
```

## Usage

### Start the Server

**Option 1: Using the start script**
```bash
./start_web.sh
```

**Option 2: Direct Python**
```bash
python3 web_interface.py
```

### Access the Interface
Open your browser and navigate to:
```
http://localhost:5000
```

## API Endpoints

### GET `/api/targets/validate`
Returns list of targets with `apphot='validate'` status.

**Response:**
```json
[
  {
    "date": "20251109",
    "target_name": "TXEri",
    "filter": "V",
    "num_images": 45,
    "submitted": 0
  }
]
```

### GET `/api/lightcurve/<date>/<target_name>/<filter>`
Returns photometry data for plotting.

**Response:**
```json
{
  "target_name": "TXEri",
  "date": "20251109",
  "filter": "V",
  "mag_type": "Standard",
  "T1": [
    {
      "index": 0,
      "hjd": 2460234.567,
      "mag": 12.345,
      "err": 0.012,
      "rejected": false
    }
  ],
  "T2": [...]
}
```

### POST `/api/target/update_status`
Updates target submission status.

**Request:**
```json
{
  "date": "20251109",
  "target_name": "TXEri",
  "filter": "V",
  "status": "approved"
}
```

**Status values:**
- `usrrej` - User rejected
- `approved` - Approved for submission
- `submitted` - Submitted to AAVSO

### POST `/api/target/save_rejected_points`
Saves rejected point indices to file.

**Request:**
```json
{
  "date": "20251109",
  "target_name": "TXEri",
  "filter": "V",
  "rejected_indices": ["T1_5", "T1_12", "T2_8"]
}
```

### POST `/api/aavso/submit`
Downloads AAVSO report file.

**Request:**
```json
{
  "date": "20251109",
  "target_name": "TXEri",
  "filter": "V"
}
```

**Response:** File download (AAVSO Extended Format)

## Database Schema

The interface interacts with the `target` table:

```sql
CREATE TABLE target (
    date TEXT NOT NULL,
    target_name TEXT NOT NULL,
    filter TEXT,
    num_images INTEGER DEFAULT 0,
    proc_bias INTEGER DEFAULT 0,
    proc_flat INTEGER DEFAULT 0,
    apphot INTEGER DEFAULT 0,
    submitted INTEGER DEFAULT 0,
    PRIMARY KEY (date, target_name, filter)
);
```

### Submission Status Values
- `0` or `NULL` - Not processed
- `'usrrej'` - User rejected
- `'approved'` - Approved, ready for AAVSO
- `'submitted'` - Submitted to AAVSO

## File Structure

```
processPhotometry/
├── web_interface.py          # Flask application
├── templates/
│   └── index.html            # Web UI with Plotly
├── start_web.sh              # Quick start script
├── wombatpipeline.db         # SQLite database
├── *_phot_transformed.tbl    # Photometry tables
├── AAVSO_AID_Report_*.txt    # AAVSO reports
└── rejections_*.json         # Rejected points
```

## Keyboard Shortcuts

- **Double-click plot**: Reset zoom
- **Shift + drag**: Pan
- **Scroll**: Zoom in/out

## Troubleshooting

### No targets showing up
- Check database: `sqlite3 wombatpipeline.db "SELECT * FROM target WHERE apphot='validate'"`
- Ensure Phase 6 transformation completed successfully
- Verify `apphot` status is set to `'validate'` (not `0` or other value)

### Lightcurve not loading
- Verify photometry table exists: `ls *_phot*.tbl`
- Check table has required columns: `HJD_UTC`, `Source_AMag_T1`, `Source_AMag_T2`
- For transformed magnitudes: `V_std_T1`, `V_std_T2`

### AAVSO report not found
- Ensure AAVSO report was generated during Phase 5
- Check for files: `ls AAVSO_AID_Report_*.txt`
- Regenerate using `pipeline.create_aavso_report()`

## Development

### Running in Debug Mode
```python
app.run(debug=True, port=5000)
```

### Custom Port
```python
app.run(debug=True, port=8080)
```

### Production Deployment
For production use, deploy with a WSGI server like Gunicorn:
```bash
pip install gunicorn
gunicorn -w 4 -b 0.0.0.0:5000 web_interface:app
```

## Security Notes

⚠️ **This interface is intended for local use only.**

For production deployment:
- Add authentication
- Use HTTPS
- Implement CSRF protection
- Validate all inputs
- Use environment variables for configuration

## License

Part of the Wombat Pipeline project.
