#!/bin/bash

# Setup script to download required AstroImageJ libraries

LIB_DIR="/Users/aegiskeller/Documents/lib"

echo "Setting up AstroImageJ libraries for Night View..."
echo "This requires downloading libraries from AstroImageJ"
echo ""
echo "Library directory: $LIB_DIR"
echo ""

# Create lib directory in Documents
mkdir -p "$LIB_DIR"

echo "ERROR: AstroImageJ libraries (ij.jar, Astronomy_.jar, Jama.jar) not found!"
echo ""
echo "Please do ONE of the following:"
echo ""
echo "1. If you have AstroImageJ installed, copy these files to $LIB_DIR:"
echo "   - ij.jar"
echo "   - Astronomy_.jar"
echo "   - Jama.jar"
echo ""
echo "2. Download AstroImageJ from https://www.astro.louisville.edu/software/astroimagej/"
echo "   Then copy the JAR files from the AstroImageJ installation to $LIB_DIR"
echo ""
echo "3. Or specify the path to your AstroImageJ installation and I'll create a symlink"

exit 1
