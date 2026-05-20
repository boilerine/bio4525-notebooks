# Bio 4525: Computational Protein Structure Analysis Module

## Description

This module is a series of five Jupyter notebooks developed for **Bio 4525: Structural Bioinformatics of Proteins** at Washington University in St. Louis. Students with no prior programming experience use Python to explore protein structure at every level of biological organization — from amino acid sequence through quaternary assembly — culminating in an integrative analysis of enzyme active sites. Each notebook combines executable code with biological narrative, guided questions, and open-ended extensions so that computation serves as a tool for biological reasoning, not an end in itself. The primary example protein throughout the series is human hemoglobin (PDB: 1HHO); the capstone notebook pivots to serine proteases to illustrate enzyme mechanism and evolutionary conservation.

## Notebooks

| # | Topic | Estimated Time | Open |
|---|-------|---------------|------|
| 1 | Primary Structure — sequence retrieval, amino acid composition, motifs | ~45 min | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/boilerine/bio4525-notebooks/blob/main/notebook_01_primary_structure.ipynb) |
| 2 | Secondary Structure — α-helices, β-sheets, DSSP assignment, visualization | ~60 min | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/boilerine/bio4525-notebooks/blob/main/notebook_02_secondary_structure.ipynb) |
| 3 | Tertiary Structure — domains, disulfide bonds, salt bridges, atomic distances | ~60 min | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/boilerine/bio4525-notebooks/blob/main/notebook_03_tertiary_structure.ipynb) |
| 4 | Quaternary Structure — subunit interfaces, buried surface area, cooperative binding | ~75 min | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/boilerine/bio4525-notebooks/blob/main/notebook_04_quaternary_structure.ipynb) |
| 5 | Integration — enzyme active sites, sequence conservation, convergent evolution | ~75 min | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/boilerine/bio4525-notebooks/blob/main/notebook_05_integration.ipynb) |

**Recommended order:** Run notebooks in sequence (1 → 5). Each notebook builds on vocabulary and concepts introduced in earlier ones.

## Prerequisites

Students should have:
- Basic knowledge of protein biochemistry (amino acids, peptide bonds, the four levels of protein structure)
- Access to a Google account (for Google Colab)
- **No prior Python or programming experience required**

## Opening in Google Colab

Each notebook runs entirely in the browser via Google Colab — no software installation needed.

1. Go to [colab.research.google.com](https://colab.research.google.com)
2. Click **File → Open notebook → GitHub** (or **Upload** if you have the file locally)
3. Open the notebook and run the setup cell (Cell 2) first
4. Follow the prompt to **Runtime → Restart session**, then **Run all**

Click any **Open in Colab** button in the Notebooks table above to launch that notebook instantly in your browser.

## Library Requirements

All libraries are installed automatically by the setup cell at the top of each notebook. No manual installation is required.

| Library | Purpose | Version |
|---------|---------|---------|
| `biopython` | Downloading and parsing PDB structure files | latest |
| `nglview` | Interactive 3-D molecular visualization | **3.0.8** (pinned) |
| `pandas` | Data tables and analysis | latest |
| `matplotlib` | Charts and plots | latest |
| `numpy` | Numerical arrays and math | latest |

> **Note:** `nglview` is pinned to version `3.0.8` because later versions have a rendering issue in Google Colab. Do not upgrade this dependency.

## Course Information

**Course:** Bio 4525 — Structural Bioinformatics of Proteins  
**Institution:** Washington University in St. Louis  
**Audience:** Junior/senior Biology, Biochemistry, and Chemistry majors  
**Enrollment:** ~36 students per year
