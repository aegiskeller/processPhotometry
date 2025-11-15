# Wombat Photometry Pipeline - System Architecture

## Overview

A complete astronomical photometry data processing system with separated backend (Python) and frontend (Java).

## Architecture

```
processPhotometry/
│
├── Backend (Python)
│   ├── pipeline.py              # Main processing pipeline
│   ├── wombatpipeline.db       # SQLite database
│   ├── data/                   # Raw FITS data
│   └── .venv/                  # Python environment
│
└── Frontend (Java/JavaFX)
    ├── pom.xml / build.gradle  # Build configs
    ├── src/main/java/          # Java source
    │   └── com/wombat/photometry/
    │       ├── PhotometryApp.java       # Main app
    │       ├── model/                   # Data models
    │       ├── database/                # DB manager
    │       └── ui/                      # UI components
    └── src/main/resources/
        └── styles/
            └── dark-theme.css   # Scientific dark theme
```

## Backend (Python)

### Technologies
- **Python 3.14** with virtual environment
- **astropy** - FITS file handling, WCS, astronomical calculations
- **numpy** - Array operations, sigma-clipped median combine
- **sqlite3** - Database for tracking
- **ASTAP** - Plate solving integration

### Pipeline Phases

**Phase 1**: Directory Scanning
- Scans `data/YYYYMMDD/` directories
- Identifies Bias, Flat, Light frames
- Populates `target` and `cal` tables

**Phase 2**: Master Calibration Creation
- Creates master bias: `mbias.fits` (3-sigma clipped median)
- Creates master flats: `mflat{filter}.fits` (normalized, sigma clipped)
- Stores metadata in `cal` table

**Phase 3**: Science Frame Processing
- Calibrates science frames (bias subtract, flat divide)
- Plate solves with ASTAP (updates FITS headers with WCS)
- Records solutions in `wcs` table

**Phase 4**: Photometry Preparation
- Queries AAVSO VSX for target information
- Calculates JD, HJD, airmass from FITS headers
- Extracts FOV from WCS solution
- Queries AAVSO VSP for comparison stars
- Updates `apphot` status ('ready', 'nocomps', 'nocal')

### Database Schema

**target** table:
```sql
date TEXT, target_name TEXT, filter TEXT, num_images INTEGER,
proc_bias TEXT, proc_flat TEXT, apphot TEXT, submitted INTEGER,
PRIMARY KEY (date, target_name, filter)
```

**cal** table:
```sql
date TEXT, cal_type TEXT, filter TEXT, num_images INTEGER,
master_path TEXT, timestamp TEXT,
PRIMARY KEY (date, cal_type, filter)
```

**wcs** table:
```sql
date TEXT, target_name TEXT, filter TEXT, filename TEXT,
success INTEGER, ra_deg REAL, dec_deg REAL, solve_time REAL,
num_stars INTEGER, faintest_mag REAL, fov_deg REAL,
error_message TEXT, timestamp TEXT,
PRIMARY KEY (date, target_name, filter, filename)
```

### Running Backend

```bash
# Full pipeline (all phases)
python -c "from pipeline import WombatPipeline; WombatPipeline('data').run()"

# Phase 3 only (calibrate & solve)
python -c "from pipeline import WombatPipeline; WombatPipeline('data').run(phase3_only=True)"

# Phase 4 only (photometry prep)
python -c "from pipeline import WombatPipeline; WombatPipeline('data').run(phase4_only=True)"
```

## Frontend (Java/JavaFX)

### Technologies
- **Java 17+**
- **JavaFX 21** - UI framework
- **SQLite JDBC** - Database connectivity
- **Maven/Gradle** - Build systems

### Features

**Dark Scientific Theme**
- High contrast (#00d4aa cyan on #0d0d0d black)
- Monospace fonts for data display
- Color-coded status indicators

**Three-Tab Interface**
1. **Targets**: Science observations with processing status
2. **WCS Solutions**: Plate solving results with coordinates
3. **Calibrations**: Master bias/flat frames

**Interactive Details**
- Click target rows to see frame-by-frame WCS solutions
- RA/Dec formatted as both degrees and HMS/DMS
- Solve statistics (stars detected, time, FOV)

**Live Statistics Dashboard**
- Total targets, ready for photometry, submitted
- Successful/failed plate solves
- Calibration frame counts

### Running Frontend

```bash
cd frontend

# Option 1: Maven
mvn clean compile && mvn javafx:run

# Option 2: Gradle  
gradle run

# Option 3: Auto-detect
./run-auto.sh
```

### Installation Requirements

**Install Maven**:
```bash
brew install maven
```

**Or Gradle**:
```bash
brew install gradle
```

**Java 17+ required**:
```bash
brew install openjdk@17
```

## Data Flow

```
Raw FITS → Phase 1 → target/cal tables
            ↓
         Phase 2 → Master calibrations (mbias.fits, mflat*.fits)
            ↓
         Phase 3 → Calibrated science frames + WCS solutions → wcs table
            ↓
         Phase 4 → VSX queries + comparison stars → apphot status
            ↓
         Frontend → Real-time database visualization
```

## UI Screenshots (Dark Theme)

### Color Scheme
- **Background**: #0d0d0d (near black)
- **Primary**: #00d4aa (cyan - headers, borders)
- **Success**: #00ff88 (green - ready status)
- **Error**: #ff4444 (red - failures)
- **Warning**: #ffa500 (orange - nocal)
- **Info**: #00a3ff (blue - normal values)

### Status Indicators
- `ready` - Green, bold (photometry ready)
- `nocal` - Orange, italic (no calibration)
- `nocomps` - Red, italic (no comparison stars)
- Date values - Blue (calibration references)
- Checkmarks - Green (submitted/solved)
- X marks - Red (failed)

## Integration

Backend and frontend are **completely decoupled**:
- Backend writes to SQLite database
- Frontend reads from same database
- No direct communication needed
- Both can run independently

Run backend to process data, then launch frontend to visualize results in real-time.

## Performance

**Backend**: Processes ~100 FITS files in ~2-3 minutes (including plate solving)
**Frontend**: Sub-second database queries, instant UI updates
**Database**: SQLite single-file, portable, no server required

## Future Enhancements

- [ ] Live backend process monitoring
- [ ] AAVSO submission integration
- [ ] Light curve plotting
- [ ] Aperture photometry UI
- [ ] Export to CSV/WebObs format
- [ ] Multi-database support
- [ ] Processing queue management
