#!/usr/bin/env python3
"""
ChemAssist Runner with API Key
"""

import os
import sys
import subprocess
import getpass

def main():
    print("🔑 ChemAssist Runner")
    print("=" * 40)
    
    # Check if API key is provided as argument
    if len(sys.argv) > 1:
        api_key = sys.argv[1]
        if api_key in ["--help", "-h", "help"]:
            print("🌐 Get your free API key from: https://console.groq.com/")
            print("\nUsage:")
            print("  python3 run_with_key.py YOUR_API_KEY")
            print("  python3 run_with_key.py  # (will prompt for key)")
            return
    else:
        print("🌐 Get your free API key from: https://console.groq.com/")
        api_key = getpass.getpass("Enter your GROQ_API_KEY: ")
    
    if not api_key or api_key == "your_groq_api_key_here":
        print("❌ Please provide a valid API key")
        print("\nUsage:")
        print("  python3 run_with_key.py YOUR_API_KEY")
        print("  python3 run_with_key.py  # (will prompt for key)")
        return
    
    # Set environment variable
    os.environ["GROQ_API_KEY"] = api_key
    
    print("✅ Setting GROQ_API_KEY and starting ChemAssist...")
    print("🌐 App will be available at: http://localhost:8501")
    print("Press Ctrl+C to stop the server")
    print("-" * 40)
    
    # Run Streamlit
    try:
        # Set PATH to include streamlit binary
        env = os.environ.copy()
        env["PATH"] = "/Users/yildiray/Library/Python/3.9/bin:" + env.get("PATH", "")
        
        subprocess.run([
            sys.executable, "-m", "streamlit", "run", "app.py", 
            "--server.port", "8501"
        ], env=env)
    except KeyboardInterrupt:
        print("\n👋 ChemAssist stopped")

if __name__ == "__main__":
    main() 