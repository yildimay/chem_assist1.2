#!/usr/bin/env python3
"""
Debug LLM response format
"""

import os
import sys
from chemassist import call_llm

def test_llm_format():
    """Test LLM response format for error fixing."""
    
    # Test input
    test_input = """
%NProcShared=4
%Mem=4GB
# B3LYP/6-31G(d)

Test calculation

0 1
C 0.0 0.0 0.0
H 1.0 0.0 0.0
H 0.0 1.0 0.0
H 0.0 0.0 1.0
"""
    
    test_log = """
Error: Basis set 6-31G(d) not found.
Please check your basis set specification.
"""
    
    system_prompt = (
        "You are an expert in computational chemistry software. "
        "A user provides a broken input file plus the tail of the log/out.\n"
        "Tasks:\n"
        "1. Diagnose the root cause (<100 words).\n"
        "2. Provide a corrected input that should run.\n"
        "IMPORTANT: You must return exactly two fenced blocks in this exact format:\n"
        "```DIAGNOSIS\n<diagnosis text here>\n```\n"
        "```FIXED_INPUT\n<complete corrected input file here>\n```\n"
        "Do not include any other text outside these blocks."
    )
    
    user_msg = f"<INPUT FILE>\n{test_input.strip()}\n<LOG TAIL>\n{test_log.strip()}"
    
    print("🧪 Testing LLM response format...")
    print("=" * 50)
    
    try:
        response = call_llm(
            system_prompt=system_prompt,
            messages=[{"role": "user", "content": user_msg}],
        )
        
        print("✅ LLM Response:")
        print(response)
        print("\n" + "=" * 50)
        
        # Test parsing
        diagnosis, fixed = "", ""
        mode = None
        
        for line in response.splitlines():
            if line.startswith("```DIAGNOSIS"):
                mode = "d"
                continue
            if line.startswith("```FIXED_INPUT"):
                mode = "f"
                continue
            if line.startswith("```"):
                mode = None
                continue
            if mode == "d":
                diagnosis += line + "\n"
            elif mode == "f":
                fixed += line + "\n"
        
        print(f"📋 Parsed Diagnosis: {len(diagnosis.strip())} chars")
        print(f"🔧 Parsed Fixed Input: {len(fixed.strip())} chars")
        
        if fixed.strip():
            print("✅ Parsing successful!")
        else:
            print("❌ Parsing failed - no FIXED_INPUT block found")
            
    except Exception as e:
        print(f"❌ Error: {e}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        os.environ["GROQ_API_KEY"] = sys.argv[1]
    test_llm_format() 