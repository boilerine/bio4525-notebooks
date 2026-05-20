# Bio 4525: Structural Bioinformatics of Proteins — Notebook Module

## Project Overview
I am building a series of 5 Jupyter notebooks for an upper-level undergraduate biology course
(Bio 4525) at Washington University in St. Louis. The course has ~36 juniors and seniors per year,
mostly General Biology, Biochemistry, and Chemistry majors. Students have NO prior Python experience.

## Pedagogical Approach — CRITICAL
- Every notebook must be beginner-friendly. Assume zero Python background.
- All code cells must have detailed comments explaining WHAT the code does and WHY.
- Use a "scaffolded" approach: early notebooks have heavily commented, partially complete code.
  Later notebooks give students more independence.
- Each notebook should have a clear learning narrative — not just code exercises, but a story
  that connects the computation to the biology.
- Include markdown cells with biological context before each code section.
- Include "Think About It" reflection questions in markdown cells throughout.
- Include a "Try It Yourself" extension exercise at the end of each notebook.

## Technical Environment
- Platform: Google Colab (cloud-based, no local installation needed)
- All notebooks must run top-to-bottom without errors in a fresh Colab session.
- Include pip install cells at the top for any non-standard libraries.
- Libraries used: BioPython, pandas, matplotlib, nglview, numpy

## Protein of Focus
- Primary example protein: Hemoglobin (PDB: 1HHO or similar)
- Secondary examples: well-characterized enzymes for Notebook 5
- Always fetch structures from RCSB PDB programmatically (not local files)

## The 5 Notebooks
1. Primary Structure — sequence retrieval, amino acid frequencies, sequence motifs
2. Secondary Structure — alpha helices, beta sheets, turns; nglview visualization
3. Tertiary Structure — domains, disulfide bonds, salt bridges, atomic distances
4. Quaternary Structure — protein-protein interfaces, surface area, cooperative binding
5. Integration — enzyme active sites, conserved catalytic residues, structure-function

## Code Style Rules
- Use descriptive variable names (not x, y, z — use `residue_count`, `helix_list`, etc.)
- Every function must have a docstring explaining its purpose, inputs, and outputs.
- Break complex operations into small, named steps with comments.
- Print intermediate results so students can see what's happening at each step.
- Avoid one-liners that are clever but hard to read.

## File Naming Convention
- notebook_01_primary_structure.ipynb
- notebook_02_secondary_structure.ipynb
- notebook_03_tertiary_structure.ipynb
- notebook_04_quaternary_structure.ipynb
- notebook_05_integration.ipynb
