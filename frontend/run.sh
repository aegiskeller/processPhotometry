#!/bin/bash
# Build and run the Wombat Photometry Viewer

cd "$(dirname "$0")"

echo "Building Wombat Photometry Viewer..."
mvn clean compile

if [ $? -eq 0 ]; then
    echo "Starting application..."
    mvn javafx:run
else
    echo "Build failed!"
    exit 1
fi
