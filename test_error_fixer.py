#!/usr/bin/env python3
"""
Test the error fixer with sample data
"""

import os
import sys
from chemassist.core.dft.error_fixer import fix_input

def test_error_fixer():
    """Test the error fixer with sample data."""
    
    # Sample problematic input
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
    
    # Sample error log
    test_log = """
Error: Basis set 6-31G(d) not found.
Please check your basis set specification.
Job ended abnormally.
"""
    
    print("🧪 Testing Error Fixer...")
    print("=" * 40)
    
    try:
        result = fix_input(
            input_text=test_input,
            log_text=test_log,
            program="gaussian"
        )
        
        print("✅ Success!")
        print(f"Diagnosis: {result['diagnosis']}")
        print(f"Fixed input length: {len(result['fixed_input'])} chars")
        print("\nFixed input:")
        print(result['fixed_input'])
        
    except Exception as e:
        print(f"❌ Error: {e}")
        print(f"Error type: {type(e).__name__}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        os.environ["GROQ_API_KEY"] = sys.argv[1]
        print(f"Using API key: {sys.argv[1][:10]}...")
    else:
        print("No API key provided. Set GROQ_API_KEY environment variable.")
        sys.exit(1)
    
    test_error_fixer() 