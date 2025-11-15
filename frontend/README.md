# Wombat Photometry Pipeline - Frontend

JavaFX-based scientific data viewer for the Wombat astronomical photometry pipeline.

## Features

- **Dark Mode Scientific Interface** - High-contrast dark theme optimized for data visualization
- **Real-time Database Monitoring** - Live connection to SQLite backend
- **Three Data Views**:
  - **Targets**: Science observation targets with processing status
  - **WCS Solutions**: Plate solving results with coordinates
  - **Calibrations**: Master bias and flat frame tracking
- **Interactive Details**: Click targets to see frame-by-frame WCS solutions
- **Statistics Dashboard**: Pipeline metrics at a glance

## Requirements

- Java 17+
- Maven 3.6+
- JavaFX 21

## Quick Start

```bash
cd frontend
./run.sh
```

Or manually:

```bash
mvn clean compile
mvn javafx:run
```

## Architecture

```
frontend/
├── pom.xml                                 # Maven configuration
├── src/main/java/com/wombat/photometry/
│   ├── PhotometryApp.java                 # Main application
│   ├── model/                             # Data models
│   │   ├── Target.java
│   │   ├── WCS.java
│   │   └── Calibration.java
│   ├── database/                          # Database layer
│   │   └── DatabaseManager.java
│   └── ui/                                # UI components
│       ├── MainView.java
│       └── StatsPanel.java
└── src/main/resources/
    └── styles/
        └── dark-theme.css                 # Dark scientific theme
```

## Database Connection

The viewer connects to `../wombatpipeline.db` by default. To use a different database, modify the path in `DatabaseManager.java`.

## Color Scheme

- **Primary Accent**: #00d4aa (Cyan)
- **Background**: #0d0d0d (Near black)
- **Success**: #00ff88 (Green)
- **Error**: #ff4444 (Red)
- **Warning**: #ffa500 (Orange)
- **Info**: #00a3ff (Blue)

## Development

The application uses:
- **JavaFX** for UI framework
- **SQLite JDBC** for database connectivity
- **Property bindings** for reactive data updates
- **Custom cell renderers** for status visualization
