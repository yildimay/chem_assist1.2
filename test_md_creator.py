#!/usr/bin/env python3
"""
Test MD input creator functionality
"""

from chemassist.core.md.input_creator import MDSpec, build_mdp, build_slurm

def test_md_creator():
    """Test all MD input creator stages."""
    
    print("🧪 Testing MD Input Creator...")
    print("=" * 50)
    
    # Create a test specification
    spec = MDSpec(
        title="Test Protein MD",
        structure_file="protein.gro",
        topology_file="protein.top",
        output_prefix="test_md",
        temperature=310.0,
        pressure=1.0,
        timestep=0.002,
        n_steps=100000,  # 200 ps
        nvt_steps=20000,  # 40 ps
        npt_steps=20000,  # 40 ps
        constraints="h-bonds",
        tcoupl="V-rescale",
        pcoupl="Parrinello-Rahman"
    )
    
    print(f"📋 Test specification created:")
    print(f"   Title: {spec.title}")
    print(f"   Temperature: {spec.temperature} K")
    print(f"   Total time: {spec.n_steps * spec.timestep:.1f} ps")
    print()
    
    # Test all stages
    stages = ["em", "nvt", "npt", "md"]
    
    for stage in stages:
        print(f"🔧 Testing {stage.upper()} stage...")
        try:
            content = build_mdp(spec, stage)
            print(f"   ✅ {stage}.mdp generated ({len(content)} chars)")
        except Exception as e:
            print(f"   ❌ Error generating {stage}.mdp: {e}")
    
    # Test SLURM script
    print(f"\n🖥️ Testing SLURM script...")
    try:
        slurm_content = build_slurm(spec)
        print(f"   ✅ SLURM script generated ({len(slurm_content)} chars)")
    except Exception as e:
        print(f"   ❌ Error generating SLURM script: {e}")
    
    print("\n🎉 All tests completed!")

if __name__ == "__main__":
    test_md_creator() 