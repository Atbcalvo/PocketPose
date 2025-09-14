# PocketPose

**Template-guided ligand posing (PyMOL + RDKit), ROSIE-ready**

> PocketPose is a small CLI that generates starting ligand poses by mapping an AlphaFold (or homology) model into a ligand-bound template pocket with PyMOL, then overlaying and MMFF-minimizing your ligand with RDKit. It outputs a pocket-aligned receptor PDB and an overlaid ligand SDF that drop straight into ROSIE/Rosetta docking. Works headless, includes a fallback for “open vs. closed” pocket mismatches, and can optionally write a pocket-only complex for quick sanity checks.

---

## Features

- **Pocket superposition:** Aligns the AF model into the template’s ligand pocket using backbone atoms (N, CA, C, O).  
- **Smart fallback:** If the AF pocket is “open” and far from the template ligand, uses mapping-based CA pairs instead of raw distance.  
- **Chem overlay:** RDKit MCS-constrained align (fallback to O3A) + ligand-only MMFF minimization.  
- **ROSIE-ready outputs:** Receptor PDB (polymer only) and overlaid ligand SDF.  
- **Headless & portable:** Prefers `pymol2`, falls back to classic PyMOL; works offscreen.  
- **Fast sanity option:** Pocket-only complex writer for quick visualization/stats (contacts/clashes).  
- **Knobs:** `--radius`, `--pocket-shell`, and `--no-complex` for speed/quality tradeoffs.

---

## Installation

You’ll need **Python 3.10+**, **RDKit**, and a headless-capable **PyMOL**.

### Linux (conda)

```bash
conda create -n posekit python=3.10 -y
conda activate posekit
conda install -c conda-forge rdkit pymol-open-source -y
python -m pip install PyQt5   # headless safety on some systems
```

### macOS (Apple Silicon) via Rosetta (x86_64 env)

```bash
conda create --platform osx-64 -n posekit_x64 python=3.10 -y
conda activate posekit_x64
conda install -c conda-forge -c schrodinger pymol-bundle rdkit -y
python -m pip install PyQt5
```

> Tip: On Apple Silicon, run commands with Rosetta when needed:
> ```bash
> arch -x86_64 python pose_from_template.py --help
> ```

---

## Quickstart

1) **Get a template PDB** (with bound ligand), an **AF model PDB**, and your **ligand SDF**.  
2) Run PocketPose (fast path skips complex assembly):

```bash
# Linux
QT_QPA_PLATFORM=offscreen PYMOL_GL_DISABLE=1 PYTHONNOUSERSITE=1 python pose_from_template.py   --template 3VXC.pdb   --af Athe_0181_3VXC.pdb   --ligand xylotetraose.sdf   --out demo_xtetra   --radius 14   --no-complex
```

**Outputs**
- `demo_xtetra_receptor.pdb` — AF protein pocket-aligned to the template frame  
- `demo_xtetra_ligand_overlaid.sdf` — ligand aligned/minimized in that frame  
- (optional) `demo_xtetra_complex.pdb` — pocket-only complex for sanity checks

Upload the receptor PDB + ligand SDF to **ROSIE/Rosetta** as starting pose inputs.

---

## Usage

```bash
python pose_from_template.py --template <TEMPLATE.pdb> --af <AF.pdb>   --ligand <LIGAND.sdf> --out <PREFIX> [--radius 10] [--no-complex] [--pocket-shell 20]
```

**Arguments**
- `--template` : Template PDB **with bound ligand** (e.g., `4G68.pdb`)  
- `--af`       : AlphaFold or homology model PDB for your target  
- `--ligand`   : Ligand SDF to place (3D preferred; 2D will be embedded)  
- `--out`      : Output prefix (files will be `<out>_*.{pdb,sdf}`)  
- `--radius`   : Pocket radius (Å) around template ligand for pocket align (default `10`)  
- `--no-complex` : Skip writing the complex PDB (faster)  
- `--pocket-shell` : When writing complex, keep protein within this Å of ligand (default `20`)

---

## How it works (under the hood)

1. **Pocket alignment (PyMOL):**
   - Rough Cα alignment to coarsely co-locate structures.  
   - Detect the **largest organic ligand** in the template.  
   - Build backbone pocket selections (N, CA, C, O) around the template ligand.  
   - If the AF pocket is too far (“open” state), **fallback** to sequence/structure mapping to select the AF pocket.  
   - Pocket-only `super` aligns AF → template pocket and writes `*_receptor.pdb`.

2. **Ligand overlay (RDKit):**
   - Save template ligand to SDF.  
   - Compute **MCS** between template ligand and your ligand; align using MCS atom map.  
   - Fallback to **O3A shape/chem** overlay if MCS is small.  
   - **MMFF minimize** the ligand (ligand-only).  
   - Write `*_ligand_overlaid.sdf`.

3. **(Optional) Pocket complex:**
   - Merge ligand + receptor pocket region and write `*_complex.pdb`.  
   - Report quick stats: contacts (≤4 Å) and hard clashes (<2 Å).

---

## ROSIE/Rosetta notes

- Some servers accept **separate** receptor (PDB) + ligand (SDF) — you’re done.  
- If a **single complex PDB** is required, either:
  - Run without `--no-complex`, **or**
  - In PyMOL: `load *_receptor.pdb; load *_ligand_overlaid.sdf; create complex, all; save complex.pdb, complex`.

Protonation: SDF contains hydrogens; Rosetta will handle standard protonation and packing during docking/refinement.

---

## Troubleshooting

- **“PyMOL unavailable in this interpreter” / `libxml2.dylib` errors**  
  Ensure you’re using the **conda environment that installed PyMOL**:
  ```bash
  conda activate posekit   # or posekit_x64 on macOS/ARM
  python -c "import pymol, pymol2; print('ok')"
  ```
  On macOS/ARM, prefer `pymol-bundle` in an **osx-64** env and run with `arch -x86_64`.

- **“Pocket selection too small”**  
  Use a larger `--radius` (e.g., 14–20). If still small, the AF model may be open; the script’s **mapping fallback** will try to recover. If it still fails, pick a **closer template**.

- **Ligand far from pocket**  
  Check `overlay RMS` in the log; values <1.5 Å are usually good. Very large values suggest the ligand is chemically dissimilar to the template ligand (MCS tiny). Try the O3A fallback (automatic) or a more similar template ligand.

- **Headless errors (`QOpenGLWidget` or Qt)**  
  Make sure these are set (or exported in your shell):
  ```bash
  export QT_QPA_PLATFORM=offscreen
  export PYMOL_GL_DISABLE=1
  ```

---

## Example commands

Download a PDB from RCSB and run:

```bash
curl -O https://files.rcsb.org/download/4G68.pdb
QT_QPA_PLATFORM=offscreen PYMOL_GL_DISABLE=1 PYTHONNOUSERSITE=1 python pose_from_template.py   --template 4G68.pdb   --af my_target_AF.pdb   --ligand my_ligand.sdf   --out myjob --radius 14 --no-complex
```

---

## License

MIT. See `LICENSE`.

---

## Acknowledgments

- **PyMOL** (open-source & `pymol-bundle`)  
- **RDKit** (BSD-3)  
- Inspired by standard template-based seeding workflows commonly used before docking.
