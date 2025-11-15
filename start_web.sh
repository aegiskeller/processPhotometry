#!/bin/bash
# Start the Wombat Pipeline Web Interface

cd "$(dirname "$0")"

# Activate virtual environment
source .venv/bin/activate

echo "Starting Wombat Pipeline Web Interface..."
echo "Navigate to: http://localhost:8080"
echo ""
echo "Press Ctrl+C to stop the server"
echo ""

python3 web_interface.py
