# Installation Guide

## Prerequisites

### Java 17 or later
Check if Java is installed:
```bash
java -version
```

If not installed:
```bash
# macOS
brew install openjdk@17

# Add to PATH
echo 'export PATH="/opt/homebrew/opt/openjdk@17/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc
```

### Maven (Recommended)
```bash
brew install maven
```

**OR**

### Gradle (Alternative)
```bash
brew install gradle
```

## Building and Running

### Option 1: Maven (Recommended)
```bash
cd frontend
mvn clean compile
mvn javafx:run
```

### Option 2: Gradle
```bash
cd frontend
gradle run
```

### Option 3: Auto-detect
```bash
cd frontend
./run-auto.sh
```

## Troubleshooting

### "command not found: mvn"
Maven is not installed. Install it:
```bash
brew install maven
```

### JavaFX errors
Make sure you have Java 17+ with JavaFX support. On macOS with brew:
```bash
brew install openjdk@17
```

### Database connection issues
Ensure the Python backend has created `wombatpipeline.db`:
```bash
cd ..
python -c "from pipeline import WombatPipeline; WombatPipeline('data').run()"
```

### Module errors
If you see module-related errors, your Java installation may need JavaFX modules. The Maven/Gradle configs handle this automatically.

## Project Structure

```
frontend/
├── pom.xml                    # Maven build (primary)
├── build.gradle              # Gradle build (alternative)
├── run.sh                    # Maven runner
├── run-auto.sh               # Auto-detect build tool
├── src/
│   └── main/
│       ├── java/
│       │   └── com/wombat/photometry/
│       │       ├── PhotometryApp.java
│       │       ├── model/
│       │       ├── database/
│       │       └── ui/
│       └── resources/
│           └── styles/
│               └── dark-theme.css
└── README.md
```

## Features Overview

### Main Window
- **Header**: Title, statistics, and control buttons
- **Tabs**: Three views for different data types
- **Footer**: Connection status

### Targets Tab
- Lists all science observation targets
- Shows processing status (bias, flat, photometry)
- Click a row to see detailed WCS solutions

### WCS Solutions Tab
- All plate solving results
- RA/Dec coordinates in degrees
- Solve statistics (stars, time, FOV)

### Calibrations Tab
- Master bias and flat frames
- Frame counts and paths
- Timestamps

### Color Coding
- **Cyan (#00d4aa)**: Headers, primary UI elements
- **Green (#00ff88)**: Success, ready status
- **Orange (#ffa500)**: Warnings, nocal status
- **Red (#ff4444)**: Errors, failures
- **Blue (#00a3ff)**: Normal status values

## Keyboard Shortcuts

- **⌘+R**: Refresh data (when focused on button)
- **⌘+Q**: Quit application
- **↑/↓**: Navigate table rows
- **Tab**: Switch between UI elements
