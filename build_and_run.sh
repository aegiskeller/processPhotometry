#!/bin/bash

# Build and run script for NightView photometry tool

# Set paths - use absolute paths to avoid directory issues
DOCS_LIB="/Users/aegiskeller/Documents/lib"
LOCAL_JAR="build/jar/nightview.jar"

# Check if JAR exists and libraries are available
if [ ! -f "$LOCAL_JAR" ]; then
    echo "Error: nightview.jar not found at $LOCAL_JAR"
    echo "You need to build the project first"
    exit 1
fi

if [ ! -d "$DOCS_LIB" ]; then
    echo "Error: Library directory not found at $DOCS_LIB"
    echo "Run setup_libs.sh first to create symlinks to AstroImageJ libraries"
    exit 1
fi

# Run the application
echo "Starting Night View application..."
echo "Classpath: $LOCAL_JAR:$DOCS_LIB/ij.jar:$DOCS_LIB/Astronomy_.jar:$DOCS_LIB/Jama.jar"
echo ""

java -cp "$LOCAL_JAR:$DOCS_LIB/ij.jar:$DOCS_LIB/Astronomy_.jar:$DOCS_LIB/Jama.jar" nightview.NightViewApp
