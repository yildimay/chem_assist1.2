#!/bin/bash

# ChemAssist Runner with API Key
echo "🔑 ChemAssist Runner"
echo "=================="

# Check if API key is provided as argument
if [ -z "$1" ]; then
    echo "❌ Please provide your API key as an argument"
    echo ""
    echo "Usage:"
    echo "  ./run_with_key.sh YOUR_GROQ_API_KEY"
    echo ""
    echo "Example:"
    echo "  ./run_with_key.sh gsk_abc123..."
    echo ""
    echo "Get your free API key from: https://console.groq.com/"
    exit 1
fi

# Set the API key and run Streamlit
echo "✅ Setting GROQ_API_KEY and starting ChemAssist..."
echo "🌐 App will be available at: http://localhost:8501"
echo ""

GROQ_API_KEY="$1" python3 -m streamlit run app.py --server.port 8501 