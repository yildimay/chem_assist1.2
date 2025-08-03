from __future__ import annotations

"""MD • Input Creator UI

Form-based generator for GROMACS .mdp files and SLURM scripts.
Supports energy minimization, NVT/NPT equilibration, and production MD.
"""

import streamlit as st
from chemassist.core.md.input_creator import MDSpec, build_mdp, build_slurm
from chemassist.utils.file_io import create_md_zip_archive

# ───────────────────────────────────────────────────────────────
# Streamlit page
# ───────────────────────────────────────────────────────────────

def show_page() -> None:  # noqa: D401
    st.header("🧮 MD • Input Creator")
    st.write("Generate GROMACS MD input files and SLURM scripts for molecular dynamics simulations.")

    with st.expander("📖 Instructions", expanded=False):
        st.markdown(
            """
            **This tool generates:**
            - Energy minimization (.mdp)
            - NVT equilibration (.mdp) 
            - NPT equilibration (.mdp)
            - Production MD (.mdp)
            - SLURM submission script
            
            **Required files:**
            - Structure file (.gro, .pdb)
            - Topology file (.top)
            
            **Workflow:**
            1. Fill out the form below
            2. Download the generated files
            3. Submit with: `sbatch run_md.sh`
            """
        )

    # ───────────────────────────────────────────────────────────────
    # Form sections
    # ───────────────────────────────────────────────────────────────

    with st.form("md_input_form"):
        st.subheader("📋 Simulation Setup")
        
        col1, col2 = st.columns(2)
        with col1:
            title = st.text_input("Simulation Title", value="MD Simulation")
            structure_file = st.text_input("Structure File", value="protein.gro")
            topology_file = st.text_input("Topology File", value="protein.top")
        
        with col2:
            output_prefix = st.text_input("Output Prefix", value="md")
            temperature = st.number_input("Temperature (K)", value=300.0, min_value=0.0, max_value=1000.0)
            pressure = st.number_input("Pressure (bar)", value=1.0, min_value=0.1, max_value=100.0)

        st.subheader("⏱️ Simulation Parameters")
        
        col1, col2, col3 = st.columns(3)
        with col1:
            timestep = st.number_input("Time Step (ps)", value=0.002, min_value=0.001, max_value=0.01, step=0.001)
            n_steps = st.number_input("Production Steps", value=50000, min_value=1000, max_value=1000000)
        
        with col2:
            nvt_steps = st.number_input("NVT Steps", value=10000, min_value=1000, max_value=100000)
            npt_steps = st.number_input("NPT Steps", value=10000, min_value=1000, max_value=100000)
        
        with col3:
            constraints = st.selectbox("Constraints", ["h-bonds", "all-bonds", "none"], index=0)
            st.write(f"Total simulation time: {(n_steps * timestep):.1f} ps")

        st.subheader("📊 Output Control")
        
        col1, col2 = st.columns(2)
        with col1:
            nstxout = st.number_input("Coordinate Output", value=1000, min_value=100, max_value=10000)
            nstvout = st.number_input("Velocity Output", value=1000, min_value=100, max_value=10000)
            nstfout = st.number_input("Force Output", value=1000, min_value=100, max_value=10000)
        
        with col2:
            nstlog = st.number_input("Log Output", value=1000, min_value=100, max_value=10000)
            nstenergy = st.number_input("Energy Output", value=1000, min_value=100, max_value=10000)

        st.subheader("🌡️ Thermostat & Barostat")
        
        col1, col2 = st.columns(2)
        with col1:
            tcoupl = st.selectbox("Temperature Coupling", ["V-rescale", "Berendsen", "no"], index=0)
            tau_t = st.number_input("Temperature τ (ps)", value=0.1, min_value=0.01, max_value=10.0, step=0.01)
        
        with col2:
            pcoupl = st.selectbox("Pressure Coupling", ["Parrinello-Rahman", "Berendsen", "no"], index=0)
            tau_p = st.number_input("Pressure τ (ps)", value=2.0, min_value=0.1, max_value=10.0, step=0.1)

        st.subheader("🔧 Cutoffs")
        
        col1, col2, col3 = st.columns(3)
        with col1:
            rlist = st.number_input("Neighbor List (nm)", value=1.0, min_value=0.5, max_value=2.0, step=0.1)
        with col2:
            rvdw = st.number_input("VDW Cutoff (nm)", value=1.0, min_value=0.5, max_value=2.0, step=0.1)
        with col3:
            rcoulomb = st.number_input("Coulomb Cutoff (nm)", value=1.0, min_value=0.5, max_value=2.0, step=0.1)

        st.subheader("🖥️ SLURM Settings")
        
        col1, col2 = st.columns(2)
        with col1:
            slurm_nodes = st.number_input("Nodes", value=1, min_value=1, max_value=10)
            slurm_ntasks_per_node = st.number_input("Tasks per Node", value=4, min_value=1, max_value=32)
            slurm_mem_per_cpu = st.number_input("Memory per CPU (GB)", value=2, min_value=1, max_value=16)
        
        with col2:
            slurm_time = st.text_input("Time Limit", value="02:00:00")
            slurm_partition = st.text_input("Partition", value="compute")

        submitted = st.form_submit_button("🚀 Generate Files")

    # ───────────────────────────────────────────────────────────────
    # Generate and display files
    # ───────────────────────────────────────────────────────────────

    if submitted:
        # Create MDSpec object
        spec = MDSpec(
            title=title,
            structure_file=structure_file,
            topology_file=topology_file,
            output_prefix=output_prefix,
            temperature=temperature,
            pressure=pressure,
            timestep=timestep,
            n_steps=n_steps,
            nvt_steps=nvt_steps,
            npt_steps=npt_steps,
            nstxout=nstxout,
            nstvout=nstvout,
            nstfout=nstfout,
            nstlog=nstlog,
            nstenergy=nstenergy,
            constraints=constraints,
            tcoupl=tcoupl,
            tau_t=tau_t,
            pcoupl=pcoupl,
            tau_p=tau_p,
            rlist=rlist,
            rvdw=rvdw,
            rcoulomb=rcoulomb,
            slurm_nodes=slurm_nodes,
            slurm_ntasks_per_node=slurm_ntasks_per_node,
            slurm_mem_per_cpu=slurm_mem_per_cpu,
            slurm_time=slurm_time,
            slurm_partition=slurm_partition
        )

        st.success("✅ Files generated successfully!")

        # Display summary
        st.subheader("📊 Simulation Summary")
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Total Time", f"{n_steps * timestep:.1f} ps")
        with col2:
            st.metric("Temperature", f"{temperature} K")
        with col3:
            st.metric("Pressure", f"{pressure} bar")

        # Generate and display files
        st.subheader("📁 Generated Files")

        # Energy minimization
        with st.expander("⚡ Energy Minimization (em.mdp)", expanded=True):
            em_content = build_mdp(spec, "em")
            st.code(em_content, language="text")
            st.download_button(
                label="💾 Download em.mdp",
                data=em_content.encode(),
                file_name="em.mdp",
                mime="text/plain"
            )

        # NVT equilibration
        with st.expander("🌡️ NVT Equilibration (nvt.mdp)", expanded=True):
            nvt_content = build_mdp(spec, "nvt")
            st.code(nvt_content, language="text")
            st.download_button(
                label="💾 Download nvt.mdp",
                data=nvt_content.encode(),
                file_name="nvt.mdp",
                mime="text/plain"
            )

        # NPT equilibration
        with st.expander("📊 NPT Equilibration (npt.mdp)", expanded=True):
            npt_content = build_mdp(spec, "npt")
            st.code(npt_content, language="text")
            st.download_button(
                label="💾 Download npt.mdp",
                data=npt_content.encode(),
                file_name="npt.mdp",
                mime="text/plain"
            )

        # Production MD
        with st.expander("🎯 Production MD (md.mdp)", expanded=True):
            md_content = build_mdp(spec, "md")
            st.code(md_content, language="text")
            st.download_button(
                label="💾 Download md.mdp",
                data=md_content.encode(),
                file_name="md.mdp",
                mime="text/plain"
            )

        # SLURM script
        with st.expander("🖥️ SLURM Script (run_md.sh)", expanded=True):
            slurm_content = build_slurm(spec)
            st.code(slurm_content, language="bash")
            st.download_button(
                label="💾 Download run_md.sh",
                data=slurm_content.encode(),
                file_name="run_md.sh",
                mime="text/plain"
            )

        # Download all files as zip
        st.subheader("📦 Download All Files")
        
        try:
            zip_data = create_md_zip_archive(spec, build_mdp, build_slurm)
            zip_filename = f"{output_prefix}_md_files.zip"
            
            st.success(f"✅ ZIP archive created with {len(zip_data)} bytes")
            
            col1, col2 = st.columns([1, 2])
            with col1:
                st.download_button(
                    label="📦 Download ZIP Archive",
                    data=zip_data,
                    file_name=zip_filename,
                    mime="application/zip",
                    help="Download all MD input files as a ZIP archive"
                )
            
            with col2:
                st.info(f"Contains: em.mdp, nvt.mdp, npt.mdp, md.mdp, run_md.sh, README.md")
                
        except Exception as e:
            st.error(f"❌ Failed to create ZIP archive: {e}")
            st.info("You can still download individual files above")

        # Next steps
        st.subheader("🚀 Next Steps")
        st.markdown(
            """
            1. **Download all files** using the buttons above
            2. **Place files in your simulation directory** with your structure and topology files
            3. **Submit the job**: `sbatch run_md.sh`
            4. **Monitor progress**: `squeue -u $USER`
            5. **Check output**: `tail -f md_*.out`
            """
        )
