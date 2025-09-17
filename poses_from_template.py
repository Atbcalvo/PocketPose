#!/usr/bin/env python3
"""
poses_from_template.py (Local Alignment)

Usage:
  python poses_from_template.py --template 4G68.pdb --af Athe_0174_4G68.pdb --ligand xylotetraose.sdf --out Athe_0174_4G68_xtetra
  # Faster run (skip complex step):
  python poses_from_template.py --template 3VXC.pdb --af Athe_0181_3VXC.pdb --ligand xylotetraose.sdf --out Athe_0181_3VXC_xtetra --radius 14 --no-complex

Outputs:
  <out>_receptor.pdb               (AF protein, pocket-aligned to template)
  <out>_ligand_overlaid.sdf        (ligand aligned to template ligand, minimized; keeps explicit H)
  <out>_complex.pdb                (receptor pocket + overlaid ligand)  [skipped if --no-complex]
"""

import os
import sys
import argparse
import tempfile

# --- headless & path hygiene ---
os.environ.setdefault("PYTHONNOUSERSITE", "1")
sys.path[:] = [p for p in sys.path if "/.pyenv/" not in p]
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("PYMOL_GL_DISABLE", "1")

# --- RDKit ---
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, rdMolAlign

# --- optional: Biopython for sequence-mapped fallback ---
try:
    from Bio import pairwise2
except Exception:
    pairwise2 = None

AA3_TO1 = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G',
    'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
    'THR':'T','TRP':'W','TYR':'Y','VAL':'V','MSE':'M'
}

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
        key = (a.segi or "", a.chain or "", a.resn or "", a.resi or "")
        counts[key] += 1

    if not counts:
        raise RuntimeError("Could not enumerate ligand residues in the template.")

    (segi, chain, resn, resi), _ = max(counts.items(), key=lambda kv: kv[1])
    return _sel_from_key(tpl_name, segi, chain, resn, resi)

# --------- sequence-mapped AF pocket (preferred local fallback) ----------
def _chain_order_table(cmd, obj):
    """
    Return per-chain residue order and 1-letter seq built from CA atoms.
    { chain: {'res': [(resi_int, resi_str, resn3)], 'seq': '...'} }
    """
    space = {'rows': []}
    cmd.iterate(f"{obj} and polymer.protein and name CA",
                "rows.append((chain, resi, resn))", space=space)
    tbl = {}
    for ch, resi, resn in space['rows']:
        if ch not in tbl:
            tbl[ch] = {'res': []}
        try:
            rint = int(str(resi).strip())
        except Exception:
            rint = None
        tbl[ch]['res'].append((rint, str(resi).strip(), str(resn).upper()))
    for ch in tbl:
        rows = tbl[ch]['res']
        rows.sort(key=lambda x: (x[0] is None, x[0] if x[0] is not None else x[1]))
        tbl[ch]['seq'] = "".join(AA3_TO1.get(r[2], 'X') for r in rows)
    return tbl

def _template_pocket_residues(cmd, templ_pocket_sel):
    """Return set of (chain,resi_str) in the template pocket (CA atoms)."""
    space = {'tpl': set()}
    cmd.iterate(f"{templ_pocket_sel} and name CA",
                "tpl.add((chain, resi))", space=space)
    return {(ch, str(rs).strip()) for ch, rs in space['tpl']}

def _aln_index_map(a_aln, b_aln):
    """Map ungapped indices: index_in_A -> index_in_B."""
    ai = bi = 0
    amap = {}
    for ca, cb in zip(a_aln, b_aln):
        if ca != '-' and cb != '-':
            amap[ai] = bi
        if ca != '-':
            ai += 1
        if cb != '-':
            bi += 1
    return amap

def map_pocket_by_sequence(cmd, tpl_obj, af_obj, templ_pocket_sel, bb="name N+CA+C+O", verbose=True):
    """
    Build AF pocket selection by mapping template pocket residues with per-chain seq alignment.
    Returns PyMOL selection string (or None if mapping failed).
    """
    if pairwise2 is None:
        if verbose:
            print("      - Biopython not available; skipping sequence-mapped fallback.")
        return None

    tpl_tab = _chain_order_table(cmd, tpl_obj)
    af_tab  = _chain_order_table(cmd, af_obj)
    tpl_pocket = _template_pocket_residues(cmd, templ_pocket_sel)
    if not tpl_pocket:
        return None

    sel_parts = []
    mapped_count = 0

    for ch in sorted(set(ch for ch, _ in tpl_pocket)):
        if ch not in tpl_tab or ch not in af_tab:
            continue
        tpl_rows = tpl_tab[ch]['res']
        af_rows  = af_tab[ch]['res']
        if not tpl_rows or not af_rows:
            continue

        tpl_seq = tpl_tab[ch]['seq']
        af_seq  = af_tab[ch]['seq']
        aln = pairwise2.align.globalxx(tpl_seq, af_seq, one_alignment_only=True)
        if not aln:
            continue
        a_tpl, a_af, *_ = aln[0]
        idx_map = _aln_index_map(a_tpl, a_af)

        tpl_resi_to_pos = {rstr: i for i, (_, rstr, _) in enumerate(tpl_rows)}
        for (tpl_ch, tpl_resi_str) in sorted([x for x in tpl_pocket if x[0]==ch],
                                             key=lambda x: tpl_resi_to_pos.get(x[1], 10**9)):
            pos_tpl = tpl_resi_to_pos.get(tpl_resi_str)
            if pos_tpl is None:
                continue
            pos_af = idx_map.get(pos_tpl)
            if pos_af is None or pos_af >= len(af_rows):
                continue
            _, af_resi_str, _ = af_rows[pos_af]
            sel_parts.append(f"( {af_obj} and chain {ch} and resi {af_resi_str} and {bb} )")
            mapped_count += 1

    if mapped_count >= 10 and sel_parts:
        mapped_sel = " or ".join(sel_parts)
        if verbose:
            print(f"      - sequence-mapped AF pocket residues: {mapped_count} backbone atoms")
        return mapped_sel

    if verbose:
        print("      - sequence-mapped fallback yielded too few atoms; skipping.")
    return None

# --------- CA-pair mapped AF pocket (structural local fallback) ----------
def mapped_af_pocket(cmd, tpl_name, af_name, lig_sel, bb="name N+CA+C+O", radius=12.0):
    """
    Map template pocket residues to AF residues via CA alignment pairs.
    Returns number of atoms in 'af_pocket_map' (0 if failed).
    """
    # collect template pocket residue keys near ligand in a local space dict
    space = {'tpl_keys': set()}
    cmd.iterate(f"({tpl_name} and polymer.protein and {bb}) within {radius:.1f} of {lig_sel}",
                "tpl_keys.add((segi,chain,resn,resi))", space=space)
    tpl_keys = space['tpl_keys']

    # CA alignment to get index pairs
    try:
        cmd.align(f"{af_name} and name CA", f"{tpl_name} and name CA")
    except Exception:
        cmd.super(f"{af_name} and name CA", f"{tpl_name} and name CA")

    pairs = cmd.get_raw_alignment(f"{af_name} and name CA", f"{tpl_name} and name CA")
    af_ca = {a.index:(a.segi,a.chain,a.resn,a.resi) for a in cmd.get_model(f"{af_name} and name CA").atom}
    tp_ca = {a.index:(a.segi,a.chain,a.resn,a.resi) for a in cmd.get_model(f"{tpl_name} and name CA").atom}

    af_keys = set()
    for iaf, itpl in pairs:
        tkey = tp_ca.get(itpl)
        akey = af_ca.get(iaf)
        if tkey in tpl_keys and akey:
            af_keys.add(akey)

    if not af_keys:
        return 0

    clauses = []
    for segi, chain, resn, resi in af_keys:
        parts = [af_name]
        if segi:  parts += ["and", f"segi {segi}"]
        if chain: parts += ["and", f"chain {chain}"]
        parts    += ["and", f"resn {resn}", "and", f"resi {resi}"]
        clauses.append("(" + " ".join(parts) + ")")

    cmd.select("af_pocket_map", f"({bb}) and (" + " or ".join(clauses) + ")")
    return cmd.count_atoms("af_pocket_map")

# ---------- core PyMOL step ----------
def pymol_pocket_super_and_export(template_pdb, af_pdb, out_receptor_pdb, pocket_radius=10.0, local_first=True):
    """
    Pipeline:
      - load template + AF; rough global CA align to co-locate
      - find largest organic ligand in the template
      - build template backbone pocket around the ligand
      - LOCAL-FIRST: sequence-mapped AF pocket → CA-pair mapping → distance-based pocket
      - super AF pocket → template pocket; save receptor + template ligand SDF
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

        # template pocket
        bb = "name N+CA+C+O"
        cmd.select("templ_pocket", f"(tpl and polymer.protein and {bb}) within {pocket_radius:.1f} of lig_templ")
        n_tpl = cmd.count_atoms("templ_pocket")

        # ---- LOCAL-FIRST LOGIC ----
        def distance_based_af(bb, r):
            cmd.select("af_pocket", f"(af and polymer.protein and {bb}) within {r:.1f} of lig_templ")
            n = cmd.count_atoms("af_pocket")
            if n < 10:
                for extra in (6.0, 10.0, 14.0):
                    cmd.select("af_pocket", f"(af and polymer.protein and {bb}) within {r+extra:.1f} of lig_templ")
                    n = cmd.count_atoms("af_pocket")
                    print(f"      - expanded AF pocket radius → {r+extra:.1f} Å : {n} atoms")
                    if n >= 10:
                        break
            return n

        n_af = 0
        if local_first:
            # 1) sequence mapping (if Biopython present)
            mapped_sel = map_pocket_by_sequence(cmd, "tpl", "af", "templ_pocket", bb=bb, verbose=True)
            if mapped_sel:
                cmd.select("af_pocket", mapped_sel)
                n_af = cmd.count_atoms("af_pocket")
                print(f"      - mapped AF pocket (sequence): {n_af} atoms")

            # 2) CA-pair mapping
            if n_af < 10:
                n_map = mapped_af_pocket(cmd, "tpl", "af", "lig_templ", bb=bb, radius=pocket_radius)
                if n_map >= 10:
                    cmd.select("af_pocket", "af_pocket_map")
                    n_af = n_map
                    print(f"      - mapped AF pocket (CA-pairs): {n_map} atoms")

            # 3) distance-based as last resort
            if n_af < 10:
                n_af = distance_based_af(bb, pocket_radius+2.0)
                print(f"      - distance-based AF pocket: {n_af} atoms")
        else:
            # original order: distance → sequence → CA-pairs
            n_af = distance_based_af(bb, pocket_radius+2.0)
            print(f"      - distance-based AF pocket: {n_af} atoms")
            if n_af < 10:
                mapped_sel = map_pocket_by_sequence(cmd, "tpl", "af", "templ_pocket", bb=bb, verbose=True)
                if mapped_sel:
                    cmd.select("af_pocket", mapped_sel)
                    n_af = cmd.count_atoms("af_pocket")
                    print(f"      - mapped AF pocket (sequence): {n_af} atoms")
            if n_af < 10:
                n_map = mapped_af_pocket(cmd, "tpl", "af", "lig_templ", bb=bb, radius=pocket_radius)
                if n_map >= 10:
                    cmd.select("af_pocket", "af_pocket_map")
                    n_af = n_map
                    print(f"      - mapped AF pocket (CA-pairs): {n_map} atoms")

        if n_tpl < 10 or n_af < 10:
            raise RuntimeError("Pocket selection too small after local mapping; try increasing --radius or a closer template.")

        # precise local pocket super (moves AF)
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

# ---------- RDKit I/O ----------
def load_sdf(path):
    mol = Chem.MolFromMolFile(path, removeHs=False)
    if mol is None:
        raise RuntimeError(f"Failed to load SDF: {path}")
    if mol.GetNumConformers() == 0:
        # embed if 3D coords missing
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
        AllChem.MMFFOptimizeMolecule(mol)
    return mol

def rdkit_overlay_and_minimize(template_lig_sdf, input_lig_sdf, out_ligand_sdf, mcs_timeout=20):
    """
    Try MCS-constrained overlay; if small/fails, fall back to O3A shape/chem overlay.
    Then MMFF minimize ligand; write out_ligand_sdf with explicit hydrogens (ROSIE-friendly).
    Returns (rms, n_mcs_atoms).
    """
    ref = load_sdf(template_lig_sdf)
    prb = load_sdf(input_lig_sdf)

    # ensure explicit H for both; keep coordinates
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

    # ligand-only minimize (keeps H)
    try:
        AllChem.MMFFOptimizeMolecule(prbH, maxIters=200)
    except Exception:
        pass

    # write SDF WITH Hs for ROSIE
    w = Chem.SDWriter(out_ligand_sdf)
    w.write(prbH)
    w.close()

    return (float(rms) if rms is not None else None), int(n_mcs)

# ---------- optional complex ----------
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

# ---------- CLI ----------
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
    ap.add_argument("--no-local-first", action="store_true",
                    help="Disable local-first pocket mapping; use distance-based pocket first")
    args = ap.parse_args()

    out_receptor = f"{args.out}_receptor.pdb"
    out_ligsdf   = f"{args.out}_ligand_overlaid.sdf"
    out_complex  = f"{args.out}_complex.pdb"

    print("[1/4] Pocket superposition (AF → template) & export receptor + template ligand...")
    print(f"      - local-first mapping: {'ON' if not args.no_local_first else 'OFF (distance-first)'}")
    tpl_lig_sdf = pymol_pocket_super_and_export(
        args.template, args.af, out_receptor,
        pocket_radius=args.radius,
        local_first=(not args.no_local_first)
    )
    print(f"      - receptor written: {out_receptor}")
    print(f"      - template ligand sdf: {tpl_lig_sdf}")

    print("[2/4] RDKit overlay (MCS/O3A) and ligand MMFF minimize (keeping H)...")
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
        clashes = 0

    # clean temp file
    try:
        os.remove(tpl_lig_sdf)
    except Exception:
        pass

    print("[4/4] Done.")
    print("Files:")
    print(f"  Receptor PDB: {out_receptor}")
    print(f"  Ligand SDF  : {out_ligsdf}  (explicit H kept)")
    if not args.no_complex:
        print(f"  Complex PDB : {out_complex}")
    if clashes > 0:
        print("Note: Close contacts detected (<2.0 Å). Docking/minimization (e.g., ROSIE) will usually resolve these.")

if __name__ == "__main__":
    main()
