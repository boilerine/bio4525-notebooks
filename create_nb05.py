"""
Generator script for notebook_05_integration.ipynb
Run: python create_nb05.py
"""
import json

def code_cell(source, cell_id):
    return {
        "cell_type": "code",
        "execution_count": None,
        "id": cell_id,
        "metadata": {},
        "outputs": [],
        "source": source
    }

def markdown_cell(source, cell_id):
    return {
        "cell_type": "markdown",
        "id": cell_id,
        "metadata": {},
        "source": source
    }

cells = []

# ── nb05-cell-01 ── Title / Intro markdown ───────────────────────────────────
cells.append(markdown_cell(
"""# Notebook 5 — Integration: Enzyme Active Sites and Structure–Function Relationships

## Welcome to the Capstone Notebook

In the first four notebooks we built up a complete picture of protein structure:

| Notebook | Level | What We Studied |
|----------|-------|-----------------|
| 1 | Primary | Amino acid sequence and composition |
| 2 | Secondary | α-helices, β-sheets, and coils |
| 3 | Tertiary | 3-D folding, domains, disulfide bonds |
| 4 | Quaternary | Subunit assembly and cooperative binding |

Now we put it all together by studying **enzyme active sites** — the precise three-dimensional arrangement of a small number of amino acids that gives an enzyme its catalytic power.

### Our Model System: Serine Proteases

A **protease** is an enzyme that cuts peptide bonds to degrade or process proteins.
A **serine protease** uses a **serine residue as the nucleophilic "chemical fist"** that attacks the peptide bond.

Serine proteases include:
- **Trypsin** — digestive enzyme; cuts after Arg and Lys
- **Chymotrypsin** — digestive enzyme; cuts after large hydrophobic residues (Phe, Trp, Tyr)
- **Thrombin** — blood-clotting cascade
- **Elastase** — lung and tissue remodeling

Although these enzymes have different substrate preferences, they all share the same catalytic mechanism — a remarkable example of **conserved molecular logic**.

### The Catalytic Triad

At the heart of every serine protease is a trio of residues called the **catalytic triad**:

| Residue | Role | Why It Matters |
|---------|------|----------------|
| **Ser 195** | Nucleophile — attacks the peptide bond | Hydroxyl (-OH) becomes reactive when deprotonated |
| **His 57** | Proton shuttle | Accepts H from Ser 195; donates H to leaving group |
| **Asp 102** | Electrostatic stabilizer | Negative charge orients and stabilizes His 57 |

These residue numbers follow **chymotrypsin numbering**, used by convention for the entire family.

### Learning Objectives

By the end of this notebook you will be able to:
1. Load and explore a protein structure programmatically (review)
2. Locate specific catalytic residues by number and measure distances between them
3. Compare sequences from two related enzymes and visualize conservation
4. Run a master analysis pipeline on an unknown enzyme
5. Explain the concept of **convergent evolution** using structural evidence""",
    "nb05-cell-01"
))

# ── nb05-cell-02 ── Setup code ────────────────────────────────────────────────
cells.append(code_cell(
"""# ─────────────────────────────────────────────────────────────────────────────
# CELL 2 — Install and import everything we need
# Run this cell first, then use Runtime > Restart Runtime, then run all cells.
# ─────────────────────────────────────────────────────────────────────────────

# Install libraries not included with Google Colab by default
# The -q flag suppresses verbose output so we can read the results more easily
!pip install -q biopython nglview

# ── Standard library imports ─────────────────────────────────────────────────
import os                          # File path operations
import math                        # Square root for distance calculations
import warnings                    # Suppress minor BioPython warnings
warnings.filterwarnings('ignore')

# ── BioPython imports ────────────────────────────────────────────────────────
from Bio.PDB import PDBList, PDBParser, PPBuilder   # Structure I/O
from Bio import pairwise2                           # Sequence alignment

# ── Data and plotting imports ────────────────────────────────────────────────
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── 3-D visualization ────────────────────────────────────────────────────────
import nglview as nv
from google.colab import output
output.enable_custom_widget_manager()   # Required for nglview in Colab

print("All libraries loaded successfully!")

# ─────────────────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS — used throughout the notebook
# ─────────────────────────────────────────────────────────────────────────────

def calculate_distance(atom1, atom2):
    '''Calculate the straight-line (Euclidean) distance between two atoms.

    Parameters
    ----------
    atom1, atom2 : Bio.PDB.Atom objects

    Returns
    -------
    float : distance in Angstroms
    '''
    # Each atom stores its x, y, z coordinates as a numpy array in .coord
    diff_vector = atom1.coord - atom2.coord
    distance = math.sqrt(sum(diff_vector ** 2))
    return round(distance, 2)


def classify_residue(residue_name):
    '''Classify an amino acid residue by its side-chain chemistry.

    Parameters
    ----------
    residue_name : str  — three-letter amino acid code (e.g. 'ALA', 'GLU')

    Returns
    -------
    str : one of 'Nonpolar', 'Polar', 'Positive', 'Negative'
    '''
    nonpolar  = {'ALA', 'VAL', 'LEU', 'ILE', 'PRO', 'PHE', 'TRP', 'MET', 'GLY'}
    polar     = {'SER', 'THR', 'CYS', 'TYR', 'ASN', 'GLN', 'HIS'}
    positive  = {'LYS', 'ARG'}
    negative  = {'ASP', 'GLU'}

    if residue_name in nonpolar:
        return 'Nonpolar'
    elif residue_name in polar:
        return 'Polar'
    elif residue_name in positive:
        return 'Positive'
    elif residue_name in negative:
        return 'Negative'
    else:
        return 'Other'   # Water, ligands, modified residues, etc.

print("Helper functions defined: calculate_distance(), classify_residue()")""",
    "nb05-cell-02"
))

# ── nb05-cell-03 ── Section 1 markdown ───────────────────────────────────────
cells.append(markdown_cell(
"""---
## Section 1 — Loading and Exploring Trypsin (PDB: 1TGN)

### Why 1TGN?

**1TGN** is the structure of **trypsin inhibited by BPTI** (bovine pancreatic trypsin inhibitor), solved to 1.9 Å resolution. We use it because:
- It contains a single, complete trypsin chain with all catalytic residues present
- The resolution is high enough to precisely locate individual atoms
- It is a classic, well-characterized structure used in hundreds of research papers

### What Is "Resolution" in Crystallography?

When scientists solve a protein structure by X-ray crystallography, the **resolution** (measured in Ångströms, Å) tells us how much detail we can see:

| Resolution | What It Means |
|-----------|---------------|
| < 1.5 Å | Atomic resolution — individual atoms clearly visible |
| 1.5–2.5 Å | High resolution — side chains well-defined |
| 2.5–3.5 Å | Medium resolution — backbone clear, side chains approximate |
| > 3.5 Å | Low resolution — rough shape only |

1 Å = 0.1 nanometers = 0.0000000001 meters — about the radius of a hydrogen atom!

### Modular Code Design

In this notebook we use a **modular** approach: we define reusable functions rather than repeating code.
This is how real bioinformatics software is written.

> **Think of functions like pipette tips** — same tool, different samples.""",
    "nb05-cell-03"
))

# ── nb05-cell-04 ── load_structure function ───────────────────────────────────
cells.append(code_cell(
"""# ─────────────────────────────────────────────────────────────────────────────
# CELL 4 — Define a reusable function to load any PDB structure
# ─────────────────────────────────────────────────────────────────────────────

def load_structure(pdb_id, file_dir='.'):
    '''Download a PDB structure and return BioPython objects for it.

    Parameters
    ----------
    pdb_id   : str  — 4-character PDB accession code (e.g. '1TGN')
    file_dir : str  — directory where the downloaded file will be saved

    Returns
    -------
    structure : Bio.PDB.Structure  — the full structure object
    model     : Bio.PDB.Model      — the first model (index 0); for X-ray
                                     structures there is usually only one model
    pdb_path  : str                — file path to the downloaded PDB file
    '''
    # Step 1: Create a PDBList object to handle downloading from the RCSB
    pdb_downloader = PDBList()

    # Step 2: Download the PDB file (if not already downloaded)
    #   obsolete=False → don't look in the obsolete archive
    #   file_type='pdb' → download the legacy PDB format (not mmCIF)
    pdb_path = pdb_downloader.retrieve_pdb_file(
        pdb_id,
        obsolete=False,
        pdir=file_dir,
        file_type='pdb'
    )

    # Step 3: Create a PDBParser and load the downloaded file
    #   QUIET=True suppresses minor warnings about non-standard residues
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_path)

    # Step 4: Extract the first (and usually only) model
    model = structure[0]

    print(f"Successfully loaded {pdb_id}")
    print(f"  File: {pdb_path}")
    print(f"  Models: {len(list(structure.get_models()))}")
    print(f"  Chains: {[chain.id for chain in model.get_chains()]}")

    return structure, model, pdb_path


print("Function load_structure() defined and ready to use.")""",
    "nb05-cell-04"
))

# ── nb05-cell-05 ── Load 1TGN and print stats ─────────────────────────────────
cells.append(code_cell(
"""# ─────────────────────────────────────────────────────────────────────────────
# CELL 5 — Load trypsin (1TGN) and explore its basic properties
# ─────────────────────────────────────────────────────────────────────────────

# Load the structure using our new function
trypsin_structure, trypsin_model, trypsin_path = load_structure('1TGN')

print()
print("=" * 55)
print("TRYPSIN STRUCTURE SUMMARY (PDB: 1TGN)")
print("=" * 55)

# ── Count chains ──────────────────────────────────────────────────────────────
chain_list = list(trypsin_model.get_chains())
print(f"Number of chains: {len(chain_list)}")

for chain in chain_list:
    # Count only standard amino acid residues (ignore water, ligands)
    amino_acids = [res for res in chain.get_residues()
                   if res.get_id()[0] == ' ']   # ' ' = standard residue flag
    print(f"  Chain {chain.id}: {len(amino_acids)} amino acid residues")

# ── Count total atoms ─────────────────────────────────────────────────────────
all_atoms = list(trypsin_model.get_atoms())
print(f"Total atoms in model: {len(all_atoms)}")

# ── Amino acid composition of chain E (the trypsin chain) ────────────────────
print()
print("Amino acid composition (chain E = trypsin):")

# Collect all standard residues from chain E
trypsin_chain = trypsin_model['E']   # Chain E is the trypsin polypeptide in 1TGN
residue_names = [res.get_resname()
                 for res in trypsin_chain.get_residues()
                 if res.get_id()[0] == ' ']

# Count each residue type using a dictionary
composition = {}
for name in residue_names:
    composition[name] = composition.get(name, 0) + 1

# Sort by count (highest first) and print top 10
sorted_composition = sorted(composition.items(), key=lambda x: x[1], reverse=True)
print(f"  Total residues: {len(residue_names)}")
print(f"  Top 10 most common:")
for residue_name, count in sorted_composition[:10]:
    percentage = (count / len(residue_names)) * 100
    print(f"    {residue_name}: {count} ({percentage:.1f}%)")""",
    "nb05-cell-05"
))

# ── nb05-cell-06 ── nglview full structure ────────────────────────────────────
cells.append(code_cell(
"""# ─────────────────────────────────────────────────────────────────────────────
# CELL 6 — Visualize the full trypsin structure in 3-D
# ─────────────────────────────────────────────────────────────────────────────

# Load the structure into nglview
view = nv.show_structure_file(trypsin_path, representations=[])

# Show trypsin (chain E) as a ribbon colored by secondary structure
# colorScheme 'sstruc' → helix=purple, sheet=yellow, coil=white (standard)
view.add_representation('cartoon',
                         selection=':E',
                         color='sstruc')

# Show the inhibitor (chain I) as sticks so we can see it binding
view.add_representation('licorice',
                         selection=':I',
                         color='#FF8C00')    # Orange for contrast

# Show calcium ion and other small molecules as spheres
view.add_representation('ball+stick',
                         selection='_CA',    # Calcium atoms
                         color='#00CED1')    # Teal

# Set initial camera view
view.camera = 'perspective'
view.parameters = {'clipDist': 0}

print("Visualization key:")
print("  Purple/yellow ribbon = trypsin (chain E), colored by secondary structure")
print("  Orange sticks        = trypsin inhibitor BPTI (chain I)")
print("  Teal sphere          = calcium ion")
print()
print("Try rotating the molecule by clicking and dragging!")
view""",
    "nb05-cell-06"
))

# ── nb05-cell-07 ── Think About It Sec1 ──────────────────────────────────────
cells.append(markdown_cell(
"""### Think About It — Section 1

1. **Trypsin cuts peptide bonds after Arg and Lys residues.** Look at the amino acid composition table you printed. How many Arg (ARG) and Lys (LYS) residues does trypsin itself have? Why doesn't trypsin digest itself?

2. **The inhibitor BPTI (chain I) binds at the active site and blocks catalysis.** Based on what you know about the catalytic mechanism, where in the active site do you predict BPTI inserts itself? What type of residue do you think it presents to Ser 195?

3. **Resolution matters.** This structure is solved at 1.9 Å resolution. What does that tell you about the confidence we have in the exact position of the catalytic residues?

4. **The structure is trimmed.** Real trypsin has 223 residues, but you may count fewer in the structure. Why might some residues be missing from a crystal structure?""",
    "nb05-cell-07"
))

# ── nb05-cell-08 ── Section 2 markdown ───────────────────────────────────────
cells.append(markdown_cell(
"""---
## Section 2 — Finding the Catalytic Triad

### The Charge-Relay Mechanism

The three residues of the catalytic triad are not just near each other — they are connected by a **network of hydrogen bonds** that dramatically amplifies the reactivity of Ser 195:

```
Asp 102 ──H-bond──► His 57 ──H-bond──► Ser 195 ──► attacks substrate
(negative)          (base)              (nucleophile)
```

Step-by-step:
1. **His 57** accepts a proton from the **Ser 195** hydroxyl (-OH), making Ser 195's oxygen a much stronger nucleophile (O⁻)
2. **Asp 102** stabilizes the positive charge that builds up on His 57 during this transfer
3. The activated Ser 195 attacks the carbonyl carbon of the peptide bond
4. A tetrahedral intermediate forms, then collapses to release the first product (the N-terminal fragment)
5. Water enters and reverses the process to release the C-terminal fragment, regenerating the free enzyme

### Standard Residue Numbering

All serine proteases use **chymotrypsin numbering** by convention:

| Residue | Role | Atom(s) Involved |
|---------|------|-----------------|
| Ser **195** | Nucleophile | **OG** (Oγ — the serine hydroxyl oxygen) |
| His **57** | Proton shuttle | **NE2** (Nε2 — the imidazole nitrogen facing Ser 195) |
| Asp **102** | Stabilizer | **OD1** / **OD2** (carboxylate oxygens) |

The key functional distance is **Ser 195 OG → His 57 NE2**, which should be about **3.0–3.5 Å** — the ideal length for a strong hydrogen bond.""",
    "nb05-cell-08"
))

# ── nb05-cell-09 ── find_catalytic_triad function ─────────────────────────────
cells.append(code_cell(
"""# ─────────────────────────────────────────────────────────────────────────────
# CELL 9 — Define a function to find the catalytic triad residues
# ─────────────────────────────────────────────────────────────────────────────

def find_catalytic_triad(model, his_num=57, asp_num=102, ser_num=195):
    '''Search all chains in a model for the catalytic triad residues.

    Serine proteases use "chymotrypsin numbering" by convention.
    This function searches every chain for residues with the given
    sequence numbers and matching residue types.

    Parameters
    ----------
    model   : Bio.PDB.Model  — the structure model to search
    his_num : int            — sequence number of the histidine  (default: 57)
    asp_num : int            — sequence number of the aspartate  (default: 102)
    ser_num : int            — sequence number of the serine     (default: 195)

    Returns
    -------
    dict with keys 'His57', 'Asp102', 'Ser195'
    Each value is a Bio.PDB.Residue object, or None if not found.
    '''
    # Start with all entries as None (not found)
    triad = {'His57': None, 'Asp102': None, 'Ser195': None}

    # Loop through every chain in the model
    for chain in model.get_chains():

        # Loop through every residue in this chain
        for residue in chain.get_residues():

            # Skip non-standard residues (water, ligands, etc.)
            if residue.get_id()[0] != ' ':
                continue

            # Get the residue sequence number and three-letter name
            res_num  = residue.get_id()[1]
            res_name = residue.get_resname()

            # Check if this is the histidine of the triad
            if res_num == his_num and res_name == 'HIS':
                triad['His57'] = residue

            # Check if this is the aspartate of the triad
            elif res_num == asp_num and res_name == 'ASP':
                triad['Asp102'] = residue

            # Check if this is the serine of the triad
            elif res_num == ser_num and res_name == 'SER':
                triad['Ser195'] = residue

    return triad


print("Function find_catalytic_triad() defined.")
print("It will search all chains and return the three active-site residues.")""",
    "nb05-cell-09"
))

# ── nb05-cell-10 ── Find triad and measure distances ─────────────────────────
cells.append(code_cell(
"""# ─────────────────────────────────────────────────────────────────────────────
# CELL 10 — Find the catalytic triad in trypsin and measure key distances
# ─────────────────────────────────────────────────────────────────────────────

# Find the triad using our function
triad = find_catalytic_triad(trypsin_model)

print("Catalytic Triad Search Results:")
print("-" * 45)

for label, residue in triad.items():
    if residue is not None:
        chain_id = residue.get_parent().id   # Get the chain this residue belongs to
        res_num  = residue.get_id()[1]
        res_name = residue.get_resname()
        print(f"  {label}: {res_name} {res_num} (chain {chain_id}) ✓")
    else:
        print(f"  {label}: NOT FOUND ✗")

print()

# ── Measure Cα–Cα distances (backbone geometry) ───────────────────────────────
print("Backbone (Cα–Cα) distances between triad residues:")
print("-" * 45)

# Helper to safely get the alpha-carbon of a residue
def get_ca(residue):
    '''Return the Cα atom of a residue, or None if missing.'''
    try:
        return residue['CA']
    except KeyError:
        return None

if all(r is not None for r in triad.values()):

    ca_his = get_ca(triad['His57'])
    ca_asp = get_ca(triad['Asp102'])
    ca_ser = get_ca(triad['Ser195'])

    if ca_his and ca_asp:
        dist_his_asp = calculate_distance(ca_his, ca_asp)
        print(f"  His57 Cα  ↔  Asp102 Cα : {dist_his_asp} Å")

    if ca_his and ca_ser:
        dist_his_ser = calculate_distance(ca_his, ca_ser)
        print(f"  His57 Cα  ↔  Ser195 Cα : {dist_his_ser} Å")

    if ca_asp and ca_ser:
        dist_asp_ser = calculate_distance(ca_asp, ca_ser)
        print(f"  Asp102 Cα ↔  Ser195 Cα : {dist_asp_ser} Å")

# ── Measure the KEY functional distance: Ser195 OG → His57 NE2 ───────────────
print()
print("Key functional distance (hydrogen bond distance):")
print("-" * 45)

try:
    ser_og  = triad['Ser195']['OG']    # Serine hydroxyl oxygen (the nucleophile)
    his_ne2 = triad['His57']['NE2']    # Histidine imidazole nitrogen facing Ser
    functional_dist = calculate_distance(ser_og, his_ne2)
    print(f"  Ser195 OG ↔ His57 NE2 : {functional_dist} Å")
    print()
    if functional_dist <= 3.5:
        print("  → Within hydrogen-bond distance (≤ 3.5 Å) ✓")
        print("  → This confirms the charge-relay network is intact.")
    else:
        print("  → Distance > 3.5 Å — the inhibitor may have perturbed the geometry.")
except KeyError as e:
    print(f"  Could not find atom {e} — check the PDB file for alternate conformations.")

# ── Also measure Asp102 OD1 → His57 ND1 (the other side of the relay) ────────
try:
    asp_od1 = triad['Asp102']['OD1']
    his_nd1 = triad['His57']['ND1']
    asp_his_dist = calculate_distance(asp_od1, his_nd1)
    print(f"  Asp102 OD1 ↔ His57 ND1: {asp_his_dist} Å")
except KeyError:
    pass""",
    "nb05-cell-10"
))

# ── nb05-cell-11 ── nglview active site ───────────────────────────────────────
cells.append(code_cell(
"""# ─────────────────────────────────────────────────────────────────────────────
# CELL 11 — Visualize the catalytic triad in 3-D
# ─────────────────────────────────────────────────────────────────────────────

# Build the NGL selection strings dynamically from the found residues
def make_ngl_selection(residue):
    '''Return an NGL selection string for a single residue, e.g. "195:E".'''
    res_num  = residue.get_id()[1]
    chain_id = residue.get_parent().id
    return f"{res_num}:{chain_id}"

# Create the view
view2 = nv.show_structure_file(trypsin_path, representations=[])

# Show the full trypsin chain as a translucent gray cartoon (context)
view2.add_representation('cartoon',
                          selection=':E',
                          color='#CCCCCC',       # Light gray
                          opacity=0.4)

# Show each triad residue as colored sticks with a larger sphere for emphasis
# Color scheme: Ser=green (nucleophile), His=blue (shuttle), Asp=red (stabilizer)
if triad['Ser195'] is not None:
    ser_sel = make_ngl_selection(triad['Ser195'])
    view2.add_representation('ball+stick',
                              selection=ser_sel,
                              color='#1A9641')   # Green for the nucleophile

if triad['His57'] is not None:
    his_sel = make_ngl_selection(triad['His57'])
    view2.add_representation('ball+stick',
                              selection=his_sel,
                              color='#2166AC')   # Blue for the proton shuttle

if triad['Asp102'] is not None:
    asp_sel = make_ngl_selection(triad['Asp102'])
    view2.add_representation('ball+stick',
                              selection=asp_sel,
                              color='#D7191C')   # Red for the stabilizer

# Also show the inhibitor so students can see how it blocks the active site
view2.add_representation('licorice',
                          selection=':I',
                          color='#FF8C00',       # Orange for inhibitor
                          opacity=0.8)

view2.camera = 'perspective'

print("Active site visualization:")
print("  Green sticks  = Ser 195 (nucleophile)")
print("  Blue sticks   = His 57 (proton shuttle)")
print("  Red sticks    = Asp 102 (stabilizer)")
print("  Orange sticks = BPTI inhibitor")
print("  Gray ribbon   = rest of trypsin (context)")
print()
print("Rotate until you can see all three triad residues.")
print("Notice how they form a nearly straight line: Asp102 → His57 → Ser195")
view2""",
    "nb05-cell-11"
))

# ── nb05-cell-12 ── Think About It Sec2 ──────────────────────────────────────
cells.append(markdown_cell(
"""### Think About It — Section 2

1. **Ser 195 OG → His 57 NE2 distance.** What distance did you measure? Is it within hydrogen-bond range (≤ 3.5 Å)? What does this tell you about the readiness of the enzyme for catalysis?

2. **The "charge relay."** Asp 102 is negatively charged at physiological pH. His 57 is normally uncharged. Ser 195 is normally a poor nucleophile. Explain in plain language how connecting these three residues transforms Ser 195 into a powerful nucleophile.

3. **Mutation experiment.** If you mutated Ser 195 → Ala (removing the hydroxyl group), what would happen to enzyme activity? What if you mutated His 57 → Ala? Would the effects be the same or different?

4. **Inhibitor location.** In the visualization, where does the BPTI inhibitor sit relative to the triad? What does its position tell you about how inhibitors can be designed to block protease activity?""",
    "nb05-cell-12"
))

# ── nb05-cell-13 ── Section 3 markdown ───────────────────────────────────────
cells.append(markdown_cell(
"""---
## Section 3 — Sequence Conservation: Trypsin vs. Chymotrypsin

### Why Compare Sequences?

Trypsin and chymotrypsin perform the same chemical reaction (cleave peptide bonds) using the same catalytic mechanism (the Ser–His–Asp triad), but they cut at **different places** in a protein:

| Enzyme | Cuts after… | Why? |
|--------|-------------|------|
| Trypsin | Arg, Lys (positive) | Its "specificity pocket" is negatively charged |
| Chymotrypsin | Phe, Trp, Tyr (large, hydrophobic) | Its pocket is large and hydrophobic |

Despite this difference, they share ~**40% sequence identity** — meaning 40% of their amino acids are the same at the same positions. This is strong evidence they **evolved from a common ancestor** (divergent evolution).

### What We Expect to Find

- **High conservation** around the catalytic triad residues (His 57, Asp 102, Ser 195)
- **Low conservation** around the "specificity pocket" (where substrate recognition occurs)
- A global alignment score that reflects their shared ancestry

### The Alignment Algorithm

We will use the **Needleman–Wunsch global alignment** algorithm (`pairwise2.align.globalms`).
It finds the optimal alignment by maximizing a score based on:
- **+2** for each matched position
- **−1** for each mismatch
- **−10** penalty to open a gap (insertion or deletion)
- **−0.5** penalty for each additional gap position

These penalties strongly discourage gaps, which is appropriate when we know the sequences have similar lengths.""",
    "nb05-cell-13"
))

# ── nb05-cell-14 ── Load chymotrypsin and align sequences ────────────────────
cells.append(code_cell(
"""# ─────────────────────────────────────────────────────────────────────────────
# CELL 14 — Load chymotrypsin (4CHA) and align with trypsin
# ─────────────────────────────────────────────────────────────────────────────

def extract_sequence(model):
    '''Extract the full amino acid sequence from a structure model.

    Uses BioPython's PPBuilder which intelligently connects polypeptide
    fragments by checking for peptide-bond geometry.

    Parameters
    ----------
    model : Bio.PDB.Model

    Returns
    -------
    str : single-letter amino acid sequence (all chains concatenated)
    '''
    builder = PPBuilder()
    sequence_parts = []

    # get_peptides() returns a list of Polypeptide objects
    # (each is a continuous stretch without chain breaks)
    for polypeptide in builder.build_peptides(model):
        sequence_parts.append(str(polypeptide.get_sequence()))

    full_sequence = ''.join(sequence_parts)
    return full_sequence


# ── Load chymotrypsin ─────────────────────────────────────────────────────────
print("Loading chymotrypsin (4CHA)...")
chymo_structure, chymo_model, chymo_path = load_structure('4CHA')

# ── Extract sequences ─────────────────────────────────────────────────────────
print()
print("Extracting sequences...")
trypsin_seq = extract_sequence(trypsin_model)
chymo_seq   = extract_sequence(chymo_model)

print(f"Trypsin  (1TGN) sequence length: {len(trypsin_seq)} residues")
print(f"Chymo    (4CHA) sequence length: {len(chymo_seq)} residues")

# ── Run global pairwise alignment ─────────────────────────────────────────────
print()
print("Running Needleman–Wunsch global alignment...")

# globalms(seq1, seq2, match_score, mismatch_score, gap_open, gap_extend)
alignments = pairwise2.align.globalms(
    trypsin_seq,   # sequence 1
    chymo_seq,     # sequence 2
    2,             # reward for a match
    -1,            # penalty for a mismatch
    -10,           # penalty to open a gap
    -0.5           # penalty for each additional gap position
)

# Take the best (highest-scoring) alignment
best_alignment = alignments[0]
aln_trypsin = best_alignment.seqA   # Aligned trypsin sequence (with gaps as '-')
aln_chymo   = best_alignment.seqB   # Aligned chymotrypsin sequence (with gaps as '-')
aln_score   = best_alignment.score

# ── Calculate identity ────────────────────────────────────────────────────────
matches   = sum(1 for a, b in zip(aln_trypsin, aln_chymo) if a == b and a != '-')
aln_len   = len(aln_trypsin)
identity  = (matches / aln_len) * 100

print(f"Alignment length: {aln_len} positions")
print(f"Identical positions: {matches}")
print(f"Sequence identity: {identity:.1f}%")
print(f"Alignment score:  {aln_score:.1f}")
print()
print(f"These enzymes are {identity:.1f}% identical — strong evidence for")
print("a common ancestral serine protease (divergent evolution).")""",
    "nb05-cell-14"
))

# ── nb05-cell-15 ── Print text alignment ─────────────────────────────────────
cells.append(code_cell(
"""# ─────────────────────────────────────────────────────────────────────────────
# CELL 15 — Print a readable alignment showing conservation
# ─────────────────────────────────────────────────────────────────────────────

def print_alignment_block(seq_a, seq_b, label_a='Trypsin ', label_b='ChymoTr ',
                           block_width=60):
    '''Print a pairwise alignment in block format with a conservation line.

    Parameters
    ----------
    seq_a, seq_b : str  — aligned sequences (with '-' for gaps)
    label_a, label_b : str — labels for each sequence (pad to same length)
    block_width : int   — number of alignment columns per printed block
    '''
    print(f"Alignment: {label_a.strip()} vs {label_b.strip()}")
    print(f"Legend: | = identical,  . = different,  - = gap")
    print()

    for start in range(0, len(seq_a), block_width):
        end = start + block_width

        # Slice this block
        block_a = seq_a[start:end]
        block_b = seq_b[start:end]

        # Build conservation line: '|' for identical, '.' for different
        conservation = ''
        for a, b in zip(block_a, block_b):
            if a == b and a != '-':
                conservation += '|'
            else:
                conservation += '.'

        # Print the block
        print(f"{label_a}{block_a}")
        print(f"        {conservation}")
        print(f"{label_b}{block_b}")
        print()


# Print the alignment (it may be long — scroll to see all of it)
print_alignment_block(aln_trypsin, aln_chymo)""",
    "nb05-cell-15"
))

# ── nb05-cell-16 ── Alignment grid visualization ──────────────────────────────
cells.append(code_cell(
"""# ─────────────────────────────────────────────────────────────────────────────
# CELL 16 — Visualize sequence conservation as a colored grid
# ─────────────────────────────────────────────────────────────────────────────

def find_aln_pos(aligned_seq, ungapped_pos):
    '''Find the alignment index corresponding to position N in the ungapped seq.

    Parameters
    ----------
    aligned_seq  : str — aligned sequence with '-' gap characters
    ungapped_pos : int — 1-based position in the ORIGINAL (ungapped) sequence

    Returns
    -------
    int : 0-based index in the aligned sequence, or -1 if not found
    '''
    # Walk through the alignment counting non-gap characters
    count = 0
    for idx, char in enumerate(aligned_seq):
        if char != '-':
            count += 1
            if count == ungapped_pos:
                return idx
    return -1   # Position not found (beyond end of sequence)


# ── Find alignment positions of the three catalytic triad residues ────────────
# These are positions in the TRYPSIN sequence (1-based)
# We need to find which alignment column each one lands in

# Trypsin catalytic triad positions in 1TGN (chymotrypsin numbering)
# Note: PPBuilder may produce a slightly shorter sequence if some residues
# are disordered — we search a ± window to be safe
triad_positions = {'His57': 57, 'Asp102': 102, 'Ser195': 195}
triad_aln_positions = {}

for label, pos in triad_positions.items():
    aln_pos = find_aln_pos(aln_trypsin, pos)
    if aln_pos == -1:
        # Try nearby positions if exact match not found (disordered residues)
        for offset in range(1, 10):
            aln_pos = find_aln_pos(aln_trypsin, pos + offset)
            if aln_pos != -1:
                break
    triad_aln_positions[label] = aln_pos

print("Alignment positions of catalytic triad residues:")
for label, pos in triad_aln_positions.items():
    if pos != -1:
        try_aa  = aln_trypsin[pos] if pos < len(aln_trypsin) else '?'
        chym_aa = aln_chymo[pos]   if pos < len(aln_chymo)   else '?'
        conserved = 'CONSERVED ✓' if try_aa == chym_aa else 'DIFFERENT'
        print(f"  {label}: alignment col {pos}  Trypsin={try_aa}  Chymo={chym_aa}  → {conserved}")
    else:
        print(f"  {label}: not found in alignment window")

# ── Draw alignment grid (60 columns centered on His57) ───────────────────────
his_aln_pos = triad_aln_positions.get('His57', 60)
window_start = max(0, his_aln_pos - 30)
window_end   = min(len(aln_trypsin), window_start + 60)

block_try  = aln_trypsin[window_start:window_end]
block_chym = aln_chymo[window_start:window_end]
n_cols     = len(block_try)

fig, ax = plt.subplots(figsize=(16, 2.5))
ax.set_xlim(0, n_cols)
ax.set_ylim(0, 2)
ax.set_aspect('equal')

for i, (a, b) in enumerate(zip(block_try, block_chym)):
    # Color: blue if identical, white if different
    if a == b and a != '-':
        bg_color = '#AEC6E8'   # Light blue = identical
    else:
        bg_color = '#FFFFFF'   # White = different

    # Draw the box
    rect = mpatches.Rectangle((i, 0), 1, 2, facecolor=bg_color,
                                edgecolor='#AAAAAA', linewidth=0.3)
    ax.add_patch(rect)

    # Write the amino acid letters (trypsin on top, chymo on bottom)
    ax.text(i + 0.5, 1.5, a, ha='center', va='center', fontsize=6)
    ax.text(i + 0.5, 0.5, b, ha='center', va='center', fontsize=6)

# Highlight catalytic triad positions with a red outline
triad_colors = {'His57': '#D7191C', 'Asp102': '#D7191C', 'Ser195': '#1A9641'}
for label, aln_pos in triad_aln_positions.items():
    col = aln_pos - window_start
    if 0 <= col < n_cols:
        color = triad_colors.get(label, 'red')
        rect = mpatches.Rectangle((col, 0), 1, 2, facecolor='none',
                                    edgecolor=color, linewidth=2.5)
        ax.add_patch(rect)
        ax.text(col + 0.5, 2.15, label, ha='center', va='bottom',
                fontsize=6, color=color, fontweight='bold')

# Labels and formatting
ax.set_yticks([0.5, 1.5])
ax.set_yticklabels(['Chymotrypsin', 'Trypsin'], fontsize=9)
ax.set_xticks(range(0, n_cols, 10))
ax.set_xticklabels(range(window_start, window_end, 10), fontsize=7)
ax.set_xlabel('Alignment position', fontsize=10)
ax.set_title('Trypsin vs. Chymotrypsin — 60-column alignment window around His 57\\n'
             '(Blue = identical residues; outlined = catalytic triad)',
             fontsize=11)

plt.tight_layout()
plt.savefig('alignment_grid.png', dpi=150, bbox_inches='tight')
plt.show()
print("Figure saved as alignment_grid.png")""",
    "nb05-cell-16"
))

# ── nb05-cell-17 ── Think About It Sec3 ──────────────────────────────────────
cells.append(markdown_cell(
"""### Think About It — Section 3

1. **Sequence identity.** What percent sequence identity did you find between trypsin and chymotrypsin? Is this more or less than you expected for two enzymes that use the same mechanism?

2. **Catalytic triad conservation.** Are His 57, Asp 102, and Ser 195 conserved between the two sequences? What would happen to enzyme function if any of these were mutated?

3. **Divergent evolution.** The fact that these enzymes share ~40% identity and the same mechanism suggests they evolved from a common ancestral protease. What selective pressure might have driven the evolution of different substrate specificities (Arg/Lys vs. Phe/Trp/Tyr)?

4. **Gap penalties.** The alignment uses a gap-open penalty of −10 and gap-extension of −0.5. Why do we penalize gaps? What biological event do gaps in a protein alignment correspond to?""",
    "nb05-cell-17"
))

# ── nb05-cell-18 ── Section 4 markdown ───────────────────────────────────────
cells.append(markdown_cell(
"""---
## Section 4 — Master Analysis Pipeline

### Building Reusable Code

So far we have written individual functions for loading structures, finding the triad, and comparing sequences.
Now we combine them into a **single master function** that can analyze any serine protease in one call.

This is the heart of **modular programming**: once you have tested each component separately, you can compose them into more powerful tools.

### What the Pipeline Does

```
analyze_serine_protease(pdb_id)
    │
    ├─► load_structure()         → downloads and parses the PDB file
    ├─► basic stats              → counts chains, residues, atoms
    ├─► amino acid composition   → classifies residues by chemistry
    ├─► find_catalytic_triad()   → locates His, Asp, Ser
    ├─► distance measurements    → Cα–Cα and functional distances
    └─► extract_sequence()       → returns the sequence for later comparison
```

The function returns a **dictionary** of results — a key–value store — so you can easily access any result by name (e.g., `result['identity']`, `result['ser_his_dist']`).""",
    "nb05-cell-18"
))

# ── nb05-cell-19 ── Master function ──────────────────────────────────────────
cells.append(code_cell(
"""# ─────────────────────────────────────────────────────────────────────────────
# CELL 19 — Define the master analysis function
# ─────────────────────────────────────────────────────────────────────────────

def analyze_serine_protease(pdb_id, his_num=57, asp_num=102, ser_num=195,
                             file_dir='.'):
    '''Complete analysis pipeline for a serine protease structure.

    Loads the structure, computes basic statistics, locates the catalytic
    triad, measures key distances, and extracts the sequence.

    Parameters
    ----------
    pdb_id   : str  — 4-character PDB accession code
    his_num  : int  — residue number of the histidine  (default: 57)
    asp_num  : int  — residue number of the aspartate  (default: 102)
    ser_num  : int  — residue number of the serine     (default: 195)
    file_dir : str  — directory for downloaded PDB files

    Returns
    -------
    results  : dict  — all computed metrics
    model    : Bio.PDB.Model
    pdb_path : str
    sequence : str
    '''
    print(f"\\n{'='*55}")
    print(f"Analyzing {pdb_id}...")
    print(f"{'='*55}")

    # ── Step 1: Load structure ────────────────────────────────────────────────
    structure, model, pdb_path = load_structure(pdb_id, file_dir=file_dir)

    # ── Step 2: Count chains and residues ─────────────────────────────────────
    chain_list = list(model.get_chains())

    total_residues = 0
    for chain in chain_list:
        aa = [r for r in chain.get_residues() if r.get_id()[0] == ' ']
        total_residues += len(aa)

    total_atoms = len(list(model.get_atoms()))

    print(f"  Chains: {len(chain_list)} | Residues: {total_residues} | Atoms: {total_atoms}")

    # ── Step 3: Amino acid composition ───────────────────────────────────────
    all_res_names = []
    for chain in chain_list:
        for residue in chain.get_residues():
            if residue.get_id()[0] == ' ':
                all_res_names.append(residue.get_resname())

    class_counts = {'Nonpolar': 0, 'Polar': 0, 'Positive': 0, 'Negative': 0, 'Other': 0}
    for name in all_res_names:
        cls = classify_residue(name)
        class_counts[cls] = class_counts.get(cls, 0) + 1

    # ── Step 4: Find the catalytic triad ──────────────────────────────────────
    triad = find_catalytic_triad(model, his_num=his_num,
                                  asp_num=asp_num, ser_num=ser_num)

    triad_found = sum(1 for r in triad.values() if r is not None)
    print(f"  Catalytic triad: {triad_found}/3 residues found")

    # ── Step 5: Measure distances ─────────────────────────────────────────────
    ser_his_dist = None
    his_asp_dist = None

    if triad['Ser195'] and triad['His57']:
        try:
            ser_his_dist = calculate_distance(triad['Ser195']['OG'],
                                               triad['His57']['NE2'])
            print(f"  Ser{ser_num} OG ↔ His{his_num} NE2 : {ser_his_dist} Å")
        except KeyError:
            # Some structures may lack the OG or NE2 atom (e.g. modified residue)
            pass

    if triad['His57'] and triad['Asp102']:
        try:
            his_asp_dist = calculate_distance(triad['His57']['CA'],
                                               triad['Asp102']['CA'])
        except KeyError:
            pass

    # ── Step 6: Extract sequence ──────────────────────────────────────────────
    sequence = extract_sequence(model)
    print(f"  Sequence length: {len(sequence)} residues")

    # ── Assemble results dictionary ───────────────────────────────────────────
    results = {
        'pdb_id':         pdb_id,
        'n_chains':       len(chain_list),
        'n_residues':     total_residues,
        'n_atoms':        total_atoms,
        'composition':    class_counts,
        'triad_found':    triad_found,
        'triad':          triad,
        'ser_his_dist':   ser_his_dist,
        'his_asp_dist':   his_asp_dist,
        'seq_length':     len(sequence),
    }

    return results, model, pdb_path, sequence


print("Master function analyze_serine_protease() defined.")
print("Ready to run on any serine protease PDB structure.")""",
    "nb05-cell-19"
))

# ── nb05-cell-20 ── Run pipeline on both enzymes and compare ─────────────────
cells.append(code_cell(
"""# ─────────────────────────────────────────────────────────────────────────────
# CELL 20 — Run the master pipeline on trypsin and chymotrypsin; compare
# ─────────────────────────────────────────────────────────────────────────────

# Run the pipeline on both enzymes
# Note: we already loaded these structures earlier, but the function
# will just re-use the already-downloaded files
trypsin_results, _, _, trypsin_seq2 = analyze_serine_protease('1TGN')
chymo_results,   _, _, chymo_seq2   = analyze_serine_protease('4CHA')

# ── Build comparison table ────────────────────────────────────────────────────
print()
print("=" * 55)
print("COMPARISON: TRYPSIN vs. CHYMOTRYPSIN")
print("=" * 55)

comparison_data = {
    'Metric': [
        'PDB ID',
        'Total residues',
        'Total atoms',
        'Triad residues found',
        'Ser–His H-bond (Å)',
        'Sequence length',
    ],
    'Trypsin (1TGN)': [
        trypsin_results['pdb_id'],
        trypsin_results['n_residues'],
        trypsin_results['n_atoms'],
        f"{trypsin_results['triad_found']}/3",
        trypsin_results['ser_his_dist'],
        trypsin_results['seq_length'],
    ],
    'Chymotrypsin (4CHA)': [
        chymo_results['pdb_id'],
        chymo_results['n_residues'],
        chymo_results['n_atoms'],
        f"{chymo_results['triad_found']}/3",
        chymo_results['ser_his_dist'],
        chymo_results['seq_length'],
    ],
}

df_comparison = pd.DataFrame(comparison_data)
print(df_comparison.to_string(index=False))

# ── Bar chart: amino acid class composition side by side ─────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle('Amino Acid Composition: Trypsin vs. Chymotrypsin', fontsize=13)

class_order  = ['Nonpolar', 'Polar', 'Positive', 'Negative']
class_colors = ['#4393C3', '#74C476', '#D7191C', '#F4A811']

for ax, results, title in zip(
        axes,
        [trypsin_results, chymo_results],
        ['Trypsin (1TGN)', 'Chymotrypsin (4CHA)']):

    comp  = results['composition']
    total = sum(comp[k] for k in class_order)
    pcts  = [(comp[k] / total) * 100 for k in class_order]

    bars = ax.bar(class_order, pcts, color=class_colors, edgecolor='white',
                  linewidth=0.8)

    # Add percentage labels on top of each bar
    for bar, pct in zip(bars, pcts):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                f'{pct:.1f}%', ha='center', va='bottom', fontsize=9)

    ax.set_title(title, fontsize=11)
    ax.set_ylabel('Percentage of residues (%)')
    ax.set_ylim(0, max(pcts) * 1.2)
    ax.spines[['top', 'right']].set_visible(False)

plt.tight_layout()
plt.savefig('composition_comparison.png', dpi=150, bbox_inches='tight')
plt.show()
print("Figure saved as composition_comparison.png")""",
    "nb05-cell-20"
))

# ── nb05-cell-21 ── Think About It Sec4 ──────────────────────────────────────
cells.append(markdown_cell(
"""### Think About It — Section 4

1. **Composition differences.** Look at the amino acid composition bar charts. Trypsin cleaves after *positively charged* residues (Arg, Lys). Does trypsin have a higher proportion of negative residues in its binding pocket than chymotrypsin? (Hint: overall composition reflects the whole protein, not just the active site — but the trend should still be visible.)

2. **Ser–His distance comparison.** How do the Ser 195 OG → His 57 NE2 distances compare between trypsin and chymotrypsin? Are they within the expected range for a hydrogen bond?

3. **Modular design advantage.** This notebook now has a single function `analyze_serine_protease()` that can process any serine protease. What are the advantages of this design over writing the same code twice? Can you think of a case where this modular approach could fail?

4. **Applying the pipeline.** In the next section you will apply this pipeline to an unknown enzyme. Before running it: what results would you expect if the enzyme is a serine protease? What if it's *not* a serine protease?""",
    "nb05-cell-21"
))

# ── nb05-cell-22 ── Capstone exercise markdown ────────────────────────────────
cells.append(markdown_cell(
"""---
## Capstone Exercise — The Mystery Enzyme

### Your Challenge

You have been given a PDB accession code for an enzyme of unknown identity: **`1SBT`**

Your task is to use all the tools built in this notebook to answer the following questions:

1. **What protein is 1SBT?** (Load it and check the printed output)
2. **Does it have a catalytic triad?** Run `find_catalytic_triad()` with the *default* His 57 / Asp 102 / Ser 195 numbers. Does it find all three residues?
3. **If not, search for the triad at different positions.** Try `his_num=64, asp_num=32, ser_num=221` instead.
4. **Is it a serine protease?** Based on the triad, the Ser–His distance, and any other evidence you find.
5. **Is it homologous to trypsin?** Align its sequence to trypsin and compute percent identity.
6. **What does the comparison tell you about evolution?**

---

### The Twist: Convergent Evolution

> *Hint: 1SBT may not be what you expect.*

Some enzymes independently evolved the same catalytic solution — the serine protease triad — even though they share **no sequence homology** and have **completely different protein folds**.

This is called **convergent evolution**: different lineages arriving at the same functional solution independently, just as dolphins and sharks independently evolved streamlined body shapes.

| Feature | Trypsin / Chymotrypsin | Mystery Enzyme |
|---------|----------------------|----------------|
| Triad chemistry | Ser–His–Asp | ? |
| Triad positions | Ser 195, His 57, Asp 102 | ? |
| Sequence identity to trypsin | 100% / ~40% | ? |
| Evolutionary origin | Common ancestral protease | ? |

### Instructions

Use the 4 blank code cells below to:
1. Load 1SBT and print basic stats
2. Search for the catalytic triad (try both default and alternate positions)
3. Visualize the active site in nglview
4. Align to trypsin and compute percent identity

**Guiding questions to answer in a markdown cell (add one if you like):**
- What organism does 1SBT come from?
- Where are the catalytic residues located (chain and number)?
- What is the percent identity to trypsin?
- Does the similar triad geometry (Ser–His distance) suggest shared ancestry or convergent evolution?
- What is the structural/biological significance of finding the same triad in two completely unrelated proteins?""",
    "nb05-cell-22"
))

# ── nb05-cell-23 to nb05-cell-26 ── Blank student cells ───────────────────────
for i, cell_id in enumerate(["nb05-cell-23", "nb05-cell-24", "nb05-cell-25", "nb05-cell-26"]):
    cells.append(code_cell(
        f"# Student code cell {i+1} of 4\n# Add your code here\n",
        cell_id
    ))

# ── nb05-cell-27 ── Final Summary markdown ────────────────────────────────────
cells.append(markdown_cell(
"""---
## Final Summary — Five Levels of Protein Structure

Congratulations on completing all five notebooks! Here is what you have learned:

| Notebook | Level | Key Concept | Tool Used |
|----------|-------|------------|-----------|
| 1 | Primary | Sequence, composition, motifs | BioPython SeqIO, pandas |
| 2 | Secondary | Helix/sheet assignment, visualization | DSSP, nglview |
| 3 | Tertiary | Domains, disulfide bonds, atomic distances | PDB parser, matplotlib |
| 4 | Quaternary | Subunit interfaces, cooperativity | Interface analysis, Hill equation |
| 5 | Integration | Enzyme active sites, conservation, evolution | Alignment, master pipeline |

### Key Biological Themes

**Structure determines function.** The precise three-dimensional arrangement of Ser 195, His 57, and Asp 102 — positioned within 3–4 Å of each other — creates the charge-relay mechanism that makes serine proteases among the most efficient catalysts in biology.

**Conservation reveals importance.** When a residue is conserved across millions of years of evolution, it is almost always because mutations at that position are lethal. The near-perfect conservation of the catalytic triad across all serine proteases is direct evidence of this principle.

**Evolution works with what it has.** Both divergent evolution (trypsin and chymotrypsin from a common ancestor) and convergent evolution (subtilisin independently inventing the same triad) are visible in protein structures. Sequence alignment lets us distinguish between the two.

### Further Reading

| Resource | What It's For |
|----------|---------------|
| [RCSB PDB](https://www.rcsb.org) | Browse and download any protein structure |
| [UniProt](https://www.uniprot.org) | Protein sequences and functional annotation |
| [Pfam](https://www.ebi.ac.uk/interpro) | Protein domain families and alignments |
| [BLAST](https://blast.ncbi.nlm.nih.gov) | Find proteins with similar sequences |
| [PyMOL](https://pymol.org) | Advanced molecular visualization |
| BioPython Tutorial | https://biopython.org/DIST/docs/tutorial/Tutorial.html |

### Looking Ahead

The methods you have used in these notebooks — structure loading, distance measurement, sequence alignment, visualization — are the same methods used in:
- Drug discovery (finding enzyme inhibitor binding sites)
- Protein engineering (designing new enzymes)
- Evolutionary biology (reconstructing ancestral proteins)
- Clinical diagnostics (understanding disease-causing mutations)

You now have the computational foundation to explore any of these directions.

---
*Bio 4525 — Structural Bioinformatics of Proteins | Washington University in St. Louis*""",
    "nb05-cell-27"
))

# ─────────────────────────────────────────────────────────────────────────────
# Assemble and write the notebook
# ─────────────────────────────────────────────────────────────────────────────

notebook = {
    "nbformat": 4,
    "nbformat_minor": 5,
    "metadata": {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3"
        },
        "language_info": {
            "name": "python",
            "version": "3.10.0"
        },
        "colab": {
            "provenance": []
        }
    },
    "cells": cells
}

output_path = "notebook_05_integration.ipynb"
with open(output_path, "w") as f:
    json.dump(notebook, f, indent=1)

print(f"Notebook written to {output_path}")
print(f"Total cells: {len(cells)}")
