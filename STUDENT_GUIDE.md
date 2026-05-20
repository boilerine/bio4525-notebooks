# Student Guide — Bio 4525 Computational Notebooks

Welcome! This guide explains everything you need to know to use the Bio 4525 Python notebooks, even if you have never written a line of code before.

---

## What Is a Jupyter Notebook?

A Jupyter notebook is a document that mixes two things in the same file:

- **Text cells** (like this guide) — formatted explanations, tables, and questions written in plain English
- **Code cells** — blocks of Python code that you can run and see the output of immediately below

Think of it like a lab notebook where the "experiments" are lines of code, and the results — numbers, tables, graphs, 3-D molecules — appear right next to the procedure that produced them.

---

## What Is Google Colab?

Google Colab is a free service that lets you run Jupyter notebooks in your web browser without installing anything on your computer. It provides computing power (a virtual machine) in the cloud, so even basic laptops can run these notebooks.

All you need is a Google account.

---

## How to Open a Notebook

1. Go to **[colab.research.google.com](https://colab.research.google.com)**
2. Sign in with your Google account
3. Click **File → Open notebook**
4. Upload the `.ipynb` file your instructor provided, or open it from Google Drive if it has been shared with you

---

## How to Run a Cell

A **cell** is one block of code or text. To run a code cell:

- Click on the cell to select it
- Press **Shift + Enter** to run it and move to the next cell
- Or click the **▶ (play) button** that appears on the left side of the cell when you hover over it

The output (printed text, a table, a graph, or a 3-D molecule) appears directly below the cell. A number in brackets like `[3]` on the left tells you the cell has been run and was the 3rd cell executed in this session.

> **Important:** Always run cells **from top to bottom, in order.** Later cells depend on variables and functions defined in earlier ones. If you skip a cell, the ones below it may fail.

---

## The First Time You Open a Notebook

Every notebook has a **setup cell** near the top (usually Cell 2) that installs the required libraries and imports them. You must run this cell first.

After it finishes, you will see a message asking you to restart the session:

1. Go to **Runtime → Restart session** (in the menu bar at the top)
2. Then go to **Runtime → Run all** to run every cell from the beginning

You only need to do this once per session. If you close the tab and come back later, repeat these steps.

---

## What to Do If You Get an Error

Errors happen to everyone — they are a normal part of coding, not a sign that something is seriously wrong. Here is how to handle them:

**Step 1: Read the error message.**
Python error messages tell you exactly what went wrong and on which line. Look at the last line of the red box — it usually says something like `NameError: name 'x' is not defined` or `KeyError: 'A'`. This is the most useful part.

**Step 2: Check the most common causes.**
- Did you run all the cells above this one? (`NameError` often means a variable was never defined because an earlier cell was skipped.)
- Did you restart the session after the setup cell? (Library import errors are often fixed by restarting.)
- Did you accidentally edit a cell and introduce a typo?

**Step 3: Search for the error.**
Copy the last line of the error message and paste it into Google. Error messages are standardized — someone else has almost certainly encountered the same one and posted a solution.

**Step 4: Ask for help.**
If you are stuck after trying the above, write down:
- Which cell failed (the cell number)
- The full error message
- What you already tried

Then ask your instructor or a classmate.

---

## How to Save Your Work

Google Colab **does not automatically save your work permanently**. Here is how to make sure you don't lose anything:

**Option 1 — Save to Google Drive (recommended)**
- Click **File → Save a copy in Drive**
- This creates a copy in your personal Google Drive that persists after the session ends

**Option 2 — Download the notebook**
- Click **File → Download → Download .ipynb**
- Save the file to your computer

**Option 3 — Manual save**
- Press **Ctrl + S** (Windows) or **Cmd + S** (Mac) to save the current state to your Drive if you opened the notebook from Drive

> **Note:** The virtual machine Colab gives you is temporary. If you leave the tab idle for too long (about 90 minutes), the session will disconnect and all variables will be cleared. Your saved notebook file is safe, but you will need to re-run it from the top.

---

## Quick Reference

| Action | How to do it |
|--------|-------------|
| Run a cell | Click the cell, then Shift + Enter |
| Run all cells | Runtime → Run all |
| Restart the session | Runtime → Restart session |
| Stop a running cell | Click the ■ (stop) button next to the cell |
| Save to Drive | File → Save a copy in Drive |
| Add a new cell | Click + Code or + Text in the toolbar |
| Undo a change | Ctrl + Z (inside a cell) |

---

*Bio 4525 — Structural Bioinformatics of Proteins | Washington University in St. Louis*
