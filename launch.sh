#!/bin/bash

# ChemAssist Launcher
echo "⚗️ ChemAssist v2 Launcher"
echo "========================"

# Fix PATH to include streamlit
export PATH="/Users/yildiray/Library/Python/3.9/bin:$PATH"

# Check if API key is provided
if [ -n "$1" ]; then
    echo "🔑 Using provided API key"
    export GROQ_API_KEY="$1"
else
    echo "⚠️  No API key provided - AI features will be disabled"
    echo "   To enable AI features, run: ./launch.sh YOUR_API_KEY"
fi

echo "🚀 Starting ChemAssist..."
echo "🌐 App will be available at: http://localhost:8501"
echo ""

# Run the app
streamlit run app.py --server.port 8501 