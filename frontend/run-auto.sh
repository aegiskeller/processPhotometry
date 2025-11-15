#!/bin/bash
# Alternative build script using Java directly

cd "$(dirname "$0")"

echo "=== Wombat Photometry Viewer ==="
echo "Checking Java installation..."

# Check for Java
if ! command -v java &> /dev/null; then
    echo "Error: Java not found. Please install Java 17 or later."
    exit 1
fi

java -version

# Check for Maven
if command -v mvn &> /dev/null; then
    echo ""
    echo "Using Maven build..."
    mvn clean compile
    if [ $? -eq 0 ]; then
        mvn javafx:run
    else
        echo "Build failed!"
        exit 1
    fi
# Check for Gradle
elif command -v gradle &> /dev/null; then
    echo ""
    echo "Using Gradle build..."
    gradle run
else
    echo ""
    echo "Error: Neither Maven nor Gradle found."
    echo ""
    echo "Please install one of the following:"
    echo "  - Maven: brew install maven"
    echo "  - Gradle: brew install gradle"
    echo ""
    echo "Or build manually using your preferred Java build tool."
    exit 1
fi
