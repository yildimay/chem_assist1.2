#!/usr/bin/env python3
"""
Helper script to set up API keys for ChemAssist.
"""

import os
import getpass

def main():
    print("🔑 ChemAssist API Key Setup")
    print("=" * 40)
    
    print("\nChoose your LLM provider:")
    print("1. Groq (recommended - free tier available)")
    print("2. OpenAI")
    
    choice = input("\nEnter your choice (1 or 2): ").strip()
    
    if choice == "1":
        print("\n🌐 Get your free Groq API key from: https://console.groq.com/")
        api_key = getpass.getpass("Enter your GROQ_API_KEY: ")
        if api_key:
            with open(".env", "w") as f:
                f.write(f"GROQ_API_KEY={api_key}\n")
            print("✅ GROQ_API_KEY saved to .env file")
        else:
            print("❌ No API key provided")
            
    elif choice == "2":
        print("\n🌐 Get your OpenAI API key from: https://platform.openai.com/api-keys")
        api_key = getpass.getpass("Enter your OPENAI_API_KEY: ")
        if api_key:
            with open(".env", "w") as f:
                f.write(f"OPENAI_API_KEY={api_key}\n")
            print("✅ OPENAI_API_KEY saved to .env file")
        else:
            print("❌ No API key provided")
    else:
        print("❌ Invalid choice")
        return
    
    print("\n🚀 You can now run ChemAssist with:")
    print("   python3 -m streamlit run app.py")

if __name__ == "__main__":
    main() 