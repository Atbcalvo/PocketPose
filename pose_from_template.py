#!/usr/bin/env python3
"""
pose_from_template.py

Usage:
  python pose_from_template.py --template 4G68.pdb --af Athe_0174_4G68.pdb --ligand xylotetraose.sdf --out Athe_0174_4G68_xtetra
  # Faster run (skip complex step):
  python pose_from_template.py --template 3VXC.pdb --af Athe_0181_3VXC.pdb --ligand xylotetraose.sdf --out Athe_0181_3VXC_xtetra --radius 14 --no-complex

Outputs:
  <out>_receptor.pdb               (AF protein, pocket-aligned to template)
  <out>_ligand_overlaid.sdf        (ligand aligned to template ligand, minimized)
  <out>_complex.pdb                (receptor pocket + overlaid ligand)  [skipped if --no-complex]
"""

import os
import sys
import argparse
import tempfile

# --- harden imports against stray user/pyenv paths & force headless ---
os.environ.setdefault("PYTHONNOUSERSITE", "1")
sys.path[:] = [p for p in sys.path if "/.pyenv/" not in p]
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("PYMOL_GL_DISABLE", "1")

# --- RDKit ---
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, rdMolAlign

# ---------- PyMOL headless launcher (pymol2 if available, else classic) ----------
def start_pymol_headless(prefer="auto"):
    """
    Returns (pm, cmd, backend)
      - pm: object with .cmd and .stop()
      - cmd: PyMOL command API
      - backend: "pymol2" or "classic"
    """
    if prefer in ("auto", "pymol2"):
        try:
            import pymol2  # type: ignore
            pm = pymol2.PyMOL()
            pm.start(['-ckq'])  # command-line, no GUI, quiet
            return pm, pm.cmd, "pymol2"
        except Exception:
            if prefer == "pymol2":
                raise
    if prefer in ("auto", "classic"):
        try:
            from pymol import cmd as _cmd, finish_launching  # type: ignore
            finish_launching(['-ckq'])
            class _PM:
                def __init__(self, cmd): self.cmd = cmd
                def stop(self):
                    try: self.cmd.quit()
                    except Exception: pass
            return _PM(_cmd), _cmd, "classic"
        except Exception as e:
            raise SystemExit(
                "PyMOL unavailable in this interpreter. Install via conda-forge/schrodinger:\n"
                "  conda install -n posekit_x64 -c conda-forge -c schrodinger pymol-bundle\n"
                f"Original error: {e}"
            )
    raise RuntimeError("No PyMOL backend available")

# ---------- helpers ----------
def _res_key(atom):
    # group by (segi, chain, resn, resi)
    return (atom.segi or "", atom.chain or "", atom.resn or "", atom.resi or "")

def _sel_from_key(obj_name, segi, chain, resn, resi):
    parts = [obj_name]
    if segi:  parts += ["and", f"segi {segi}"]
    if chain: parts += ["and", f"chain {chain}"]
    if resn:  parts += ["and", f"resn {resn}"]
    if resi:  parts += ["and", f"resi {resi}"]
    return "(" + " ".join(parts) + ")"

def pick_template_ligand(cmd, tpl_name: str) -> str:
    """
    Return a selection string for the largest organic (non-water) ligand residue in the template object.
    """
    base_sel = f"({tpl_name} and organic and not resn HOH)"
    if cmd.count_atoms(base_sel) == 0:
        raise RuntimeError("No organic ligand found in template (is a ligand present?)")

    mdl = cmd.get_model(base_sel)
    from collections import defaultdict
    counts = defaultdict(int)
    for a in mdl.atom:
        counts[_res_key(a)] += 1

    if not counts:
        raise RuntimeError("Could not enumerate ligand residues in the template.")

    (segi, chain, resn, resi), _ = max(counts.items(), key=lambda kv: kv[1])
    return _sel_from_key(tpl_name, segi, chain, resn, resi)

def pymol_pocket_super_and_export(template_pdb, af_pdb, out_receptor_pdb, pocket_radius=10.0):
    """
    Pipeline:
      - load template + AF
      - rough global align AF→template (Cα) to get them in the same neighborhood
      - find largest organic ligand in template
      - build pocket selections (backbone) on both within radius of template ligand
      - pocket-only 'super' (moves AF precisely)
      - save AF polymer.protein as receptor PDB
      - save template ligand to temp SDF for RDKit
    """
    pm, cmd, backend = start_pymol_headless()
    try:
        cmd.reinitialize()
        cmd.load(template_pdb, "tpl")
        cmd.load(af_pdb, "af")

        # rough global align (mobile → target)
        try:
            rms_align = cmd.align("af and name CA", "tpl and name CA")[0]
        except Exception:
            rms_align = cmd.super("af and name CA", "tpl and name CA")[0]
        print(f"      - rough CA align RMS: {rms_align:.3f} Å")

        # pick template ligand
        lig_sel = pick_template_ligand(cmd, "tpl")
        cmd.select("lig_templ", lig_sel)

        # build pockets around template ligand
        bb = "name N+CA+C+O"
        cmd.select("templ_pocket", f"(tpl and polymer.protein and {bb}) within {pocket_radius:.1f} of lig_templ")
        # slightly more permissive on AF pocket
        cmd.select("af_pocket",    f"(af  and polymer.protein and {bb}) within {pocket_radius+2.0:.1f} of lig_templ")

        n_tpl = cmd.count_atoms("templ_pocket")
        n_af  = cmd.count_atoms("af_pocket")
        print(f"      - pocket sizes (tpl, af): {n_tpl}, {n_af}")

        # auto-expand AF pocket if too small
        if n_af < 10:
            for extra in (6.0, 10.0, 14.0):
                cmd.select("af_pocket", f"(af and polymer.protein and {bb}) within {pocket_radius+extra:.1f} of lig_templ")
                n_af = cmd.count_atoms("af_pocket")
                print(f"      - expanded AF pocket radius → {pocket_radius+extra:.1f} Å : {n_af} atoms")
                if n_af >= 10:
                    break

        if n_tpl < 10 or n_af < 10:
            raise RuntimeError("Pocket selection too small after rough align; try increasing --radius.")

        # precise local pocket super
        rms_pocket = cmd.super("af_pocket", "templ_pocket")[0]
        print(f"      - pocket super RMS: {rms_pocket:.3f} Å")
        cmd.delete("af_pocket")
        cmd.delete("templ_pocket")

        # save receptor and template ligand
        cmd.save(out_receptor_pdb, "af and polymer.protein")
        tmp_sdf = tempfile.NamedTemporaryFile(suffix=".sdf", delete=False).name
        cmd.save(tmp_sdf, "lig_templ")

        print(f"[PyMOL backend: {backend}]")
        return tmp_sdf
    finally:
        pm.stop()

def load_sdf(path):
    mol = Chem.MolFromMolFile(path, removeHs=False)
    if mol is None:
        raise RuntimeError(f"Failed to load SDF: {path}")
    if mol.GetNumConformers() == 0:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
        AllChem.MMFFOptimizeMolecule(mol)
    return mol

def rdkit_overlay_and_minimize(template_lig_sdf, input_lig_sdf, out_ligand_sdf, mcs_timeout=20):
    """
    Try MCS-constrained overlay; if small/fails, fall back to O3A shape/chem overlay.
    Then MMFF minimize ligand; write out_ligand_sdf.
    Returns (rms, n_mcs_atoms).
    """
    ref = load_sdf(template_lig_sdf)
    prb = load_sdf(input_lig_sdf)

    refH = Chem.AddHs(ref, addCoords=True)
    prbH = Chem.AddHs(prb, addCoords=True)

    rms = None
    n_mcs = 0
    try:
        mcs = rdFMCS.FindMCS(
            [refH, prbH],
            ringMatchesRingOnly=True,
            completeRingsOnly=True,
            timeout=mcs_timeout,
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareOrderExact,
        )
        if mcs.smartsString:
            patt = Chem.MolFromSmarts(mcs.smartsString)
            mt = refH.GetSubstructMatch(patt)
            mq = prbH.GetSubstructMatch(patt)
            if mt and mq and len(mt) == len(mq) and len(mt) >= 5:
                amap = list(zip(mq, mt))  # (probe_idx, ref_idx)
                rms = rdMolAlign.AlignMol(prbH, refH, atomMap=amap)
                n_mcs = len(amap)
    except Exception:
        pass

    if n_mcs < 5:
        try:
            ref_props = AllChem.MMFFGetMoleculeProperties(refH, mmffVariant='MMFF94')
            prb_props = AllChem.MMFFGetMoleculeProperties(prbH, mmffVariant='MMFF94')
            o3a = rdMolAlign.GetO3A(prbH, refH, prb_props=prb_props, ref_props=ref_props)
            rms = o3a.Align()
            n_mcs = 0
        except Exception:
            ref_idx = [a.GetIdx() for a in refH.GetAtoms() if a.GetAtomicNum() > 1]
            prb_idx = [a.GetIdx() for a in prbH.GetAtoms() if a.GetAtomicNum() > 1]
            common = min(len(ref_idx), len(prb_idx))
            if common > 0:
                amap = list(zip(prb_idx[:common], ref_idx[:common]))
                rms = rdMolAlign.AlignMol(prbH, refH, atomMap=amap)
            n_mcs = 0

    try:
        AllChem.MMFFOptimizeMolecule(prbH, maxIters=200)
    except Exception:
        pass

    # --- Keep explicit hydrogens for ROSIE ---
    w = Chem.SDWriter(out_ligand_sdf)
    w.write(prbH)   # prbH already has explicit H and 3D coords
    w.close()

    return (float(rms) if rms is not None else None), int(n_mcs)

def assemble_complex(receptor_pdb, ligand_sdf, out_complex_pdb, pocket_shell=20.0):
    """
    Fast pocket-only complex: save protein within pocket_shell Å of ligand + ligand.
    Also computes simple contact/clash counts against that pocket.
    """
    pm, cmd, backend = start_pymol_headless(prefer="pymol2")  # prefer pymol2 if available
    try:
        cmd.reinitialize()
        cmd.load(receptor_pdb, "receptor")
        cmd.load(ligand_sdf, "ligand")

        cmd.select("prot", "receptor and polymer.protein")
        cmd.select("prot_near", f"prot within {float(pocket_shell):.1f} of ligand")

        contacts = cmd.count_atoms("ligand within 4.0 of prot_near")
        clashes  = cmd.count_atoms("ligand within 2.0 of prot_near")

        cmd.save(out_complex_pdb, "prot_near or ligand")
        return contacts, clashes
    finally:
        pm.stop()

def main():
    ap = argparse.ArgumentParser(description="Template-guided ligand pose into AF model (PyMOL + RDKit).")
    ap.add_argument("--template", required=True, help="Template PDB with bound ligand (e.g., 4G68.pdb)")
    ap.add_argument("--af",        required=True, help="AF target PDB (e.g., Athe_0174_4G68.pdb)")
    ap.add_argument("--ligand",    required=True, help="Input ligand SDF to place (e.g., xylotetraose.sdf)")
    ap.add_argument("--out",       required=True, help="Output prefix (e.g., Athe_0174_4G68_xtetra)")
    ap.add_argument("--radius",    type=float, default=10.0, help="Pocket radius Å around template ligand (default 10)")
    ap.add_argument("--no-complex", action="store_true",
                    help="Skip complex assembly (faster; only write receptor + ligand SDF)")
    ap.add_argument("--pocket-shell", type=float, default=20.0,
                    help="If assembling complex, save only protein within this Å of ligand (default 20)")
    args = ap.parse_args()

    out_receptor = f"{args.out}_receptor.pdb"
    out_ligsdf   = f"{args.out}_ligand_overlaid.sdf"
    out_complex  = f"{args.out}_complex.pdb"

    print("[1/4] Pocket superposition (AF → template) & export receptor + template ligand...")
    tpl_lig_sdf = pymol_pocket_super_and_export(args.template, args.af, out_receptor, pocket_radius=args.radius)
    print(f"      - receptor written: {out_receptor}")
    print(f"      - template ligand sdf: {tpl_lig_sdf}")

    print("[2/4] RDKit overlay (MCS) and ligand MMFF minimize...")
    rms, n_mcs = rdkit_overlay_and_minimize(tpl_lig_sdf, args.ligand, out_ligsdf)
    if rms is not None:
        print(f"      - overlay RMS: {rms:.3f} Å (n_mcs={n_mcs})")
    else:
        print("      - overlay RMS: n/a (fallback alignment)")

    if not args.no_complex:
        print("[3/4] Assemble complex (pocket-only) for sanity...")
        contacts, clashes = assemble_complex(out_receptor, out_ligsdf, out_complex, pocket_shell=args.pocket_shell)
        print(f"      - contacts ≤4 Å: {contacts}")
        print(f"      - hard clashes <2.0 Å: {clashes}")
        print(f"      - complex written: {out_complex}")
    else:
        print("[3/4] Skipping complex assembly (--no-complex).")
        clashes = 0  # so Step 4's note logic is safe

    # clean temp file
    try:
        os.remove(tpl_lig_sdf)
    except Exception:
        pass

    print("[4/4] Done.")
    print("Files:")
    print(f"  Receptor PDB: {out_receptor}")
    print(f"  Ligand SDF  : {out_ligsdf}")
    if not args.no_complex:
        print(f"  Complex PDB : {out_complex}")
    if clashes > 0:
        print("Note: Close contacts detected (<2.0 Å). Docking/minimization (e.g., ROSIE) will usually resolve these.")

if __name__ == "__main__":
    main()
