import streamlit as st

def show_page() -> None:
    st.header(" DFT • Error Fixer")
    st.info("Upload a broken .gjf/.log duo and get a patched input file. \n\n*Implementation coming soon.*")
from chemassist.core.dft.pdf_methods import extract
...
pdf = st.file_uploader(..., type=\"pdf\")
if st.button(\"Extract methods\") and pdf:
    hits = extract(pdf.read())
    st.dataframe([h.__dict__ for h in hits])
