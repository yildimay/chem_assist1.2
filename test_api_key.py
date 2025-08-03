#!/usr/bin/env python3
"""
Test API key functionality
"""

import os
import sys
from chemassist import call_llm

def test_api_key(api_key):
    """Test if the API key works with a simple request."""
    os.environ["GROQ_API_KEY"] = api_key
    
    try:
        response = call_llm(
            messages=[{"role": "user", "content": "Say 'Hello from ChemAssist!'"}],
            model="llama3-70b-8192"
        )
        print("✅ API key works!")
        print(f"Response: {response}")
        return True
    except Exception as e:
        print(f"❌ API key failed: {e}")
        return False

def main():
    if len(sys.argv) > 1:
        api_key = sys.argv[1]
    else:
        print("Usage: python3 test_api_key.py YOUR_API_KEY")
        return
    
    print("🧪 Testing API key...")
    test_api_key(api_key)

if __name__ == "__main__":
    main() 