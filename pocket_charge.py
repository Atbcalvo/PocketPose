#!/usr/bin/env python3

"""
A tiny CLI that measures pocket volume with fpocket and computes electrostatics with PDB2PQR+APBS, then launches PyMOL headless to color the receptor surface by the APBS potential and (optionally) show the fpocket pocket surface.
Inputs: a pocket-aligned receptor PDB (and optional ligand SDF).
Outputs: an ESP-colored PyMOL session (*_esp.pse), a high-res PNG, the APBS map (*_pot.dx), and the fpocket results folder.
Open the .pse in PyMOL to explore interactively; or load the .dx and apply a blue-white-red ramp if you prefer.
"""

import os, sys, argparse, subprocess, shutil, re, json
from pathlib import Path

# Headless & clean import env
os.environ.setdefault("PYTHONNOUSERSITE", "1")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("PYMOL_GL_DISABLE", "1")

try:
    import numpy as np
except Exception:
    print("ERROR: numpy required. Install with: conda install -c conda-forge numpy", file=sys.stderr)
    sys.exit(1)

# ---------- tiny utils ----------
def run(cmd, cwd=None, check=True):
    print(">>", " ".join(cmd))
    p = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if check and p.returncode != 0:
        print(p.stdout)
        raise SystemExit(f"Command failed: {' '.join(cmd)}")
    return p.stdout

def which_or_die(name, hint):
    path = shutil.which(name)
    if not path:
        raise SystemExit(f"Required executable '{name}' not found in PATH.\n{hint}")
    return path

def load_sdf_centroid(sdf_path):
    """Very simple SDF centroid (V2000-ish)."""
    coords = []
    try:
        lines = Path(sdf_path).read_text().splitlines()
    except Exception:
        return None
    if len(lines) < 5:
        return None
    natoms = None
    try:
        natoms = int(lines[3][:3])
    except Exception:
        pass
    block = lines[4:4+natoms] if natoms else lines
    for line in block:
        toks = line.split()
        if len(toks) >= 3:
            try:
                x, y, z = float(toks[0]), float(toks[1]), float(toks[2])
                coords.append((x, y, z))
            except:
                pass
    if not coords:
        return None
    arr = np.array(coords, dtype=float)
    return arr.mean(axis=0)

def parse_pqr_xyz(pqr_file):
    xyz = []
    with open(pqr_file, "r") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                    xyz.append((x, y, z))
                except:
                    pass
    return np.array(xyz, dtype=float) if xyz else None

def pdb_bbox(pdb_path):
    xyz = []
    with open(pdb_path, "r") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                    xyz.append((x, y, z))
                except:
                    pass
    if not xyz:
        return None
    arr = np.array(xyz, dtype=float)
    return arr.min(axis=0), arr.max(axis=0)

# ---------- PyMOL helper ----------
def pymol_session_png(receptor_pdb, ligand_sdf, dx_map, pocket_vert_pqr, out_pse, out_png):
    # prefer pymol2 if available
    try:
        import pymol2
        pm = pymol2.PyMOL()
        pm.start(['-ckq'])
        cmd = pm.cmd
        backend = "pymol2"
    except Exception:
        from pymol import cmd as _cmd, finish_launching
        finish_launching(['-ckq'])
        class _PM:
            def __init__(self, cmd): self.cmd = cmd
            def stop(self):
                try: self.cmd.quit()
                except: pass
        pm = _PM(_cmd); cmd = _cmd; backend = "classic"
    print(f"[PyMOL backend: {backend}]")
    try:
        cmd.reinitialize()
        cmd.load(receptor_pdb, "rec")
        if ligand_sdf and Path(ligand_sdf).exists():
            cmd.load(ligand_sdf, "lig")
            cmd.show("sticks", "lig")
            cmd.color("tv_green", "lig")
        if pocket_vert_pqr and Path(pocket_vert_pqr).exists():
            cmd.load(pocket_vert_pqr, "pock")
            cmd.show("surface", "pock")
            cmd.set("transparency", 0.3, "pock")
            cmd.color("yellow", "pock")
        # surface & electrostatics
        cmd.show_as("surface", "rec")
        cmd.set("surface_quality", 1)
        cmd.color("white", "rec")
        if dx_map and Path(dx_map).exists():
            cmd.load(dx_map, "pot")
            # explicit red (negative) → white → blue (positive)
            cmd.ramp_new("esp_ramp", "pot", [-5, 0, +5], ["red", "white", "blue"])
            cmd.set("surface_color", "esp_ramp", "rec")
            # optional contours
            cmd.isosurface("esp_pos", "pot", +3.0); cmd.color("blue", "esp_pos")
            cmd.isosurface("esp_neg", "pot", -3.0); cmd.color("red",  "esp_neg")
        cmd.zoom("lig" if "lig" in cmd.get_object_list() else "rec", buffer=8.0)
        if out_pse:
            cmd.save(out_pse)
        try:
            cmd.bg_color("white")
            cmd.set("ray_opaque_background", 0)
            cmd.ray(2000, 1500)
        except Exception:
            pass
        if out_png:
            cmd.png(out_png, dpi=300, ray=0)
    finally:
        try:
            pm.stop()
        except Exception:
            pass

# ---------- main workflow ----------
def main():
    ap = argparse.ArgumentParser(
        description="Run fpocket (pocket volume) + APBS (electrostatics) and produce a PyMOL session/PNG with ESP-colored surface."
    )
    ap.add_argument("--receptor", required=True, help="Pocket-aligned receptor PDB (e.g., *_receptor.pdb)")
    ap.add_argument("--ligand", help="Ligand SDF near the pocket (optional but recommended)")
    ap.add_argument("--out", required=True, help="Output prefix (e.g., sample1)")
    ap.add_argument("--pH", type=float, default=7.0, help="pH for PDB2PQR (default 7.0)")
    ap.add_argument("--ff", default="amber", help="Force field for PDB2PQR (amber, CHARMM, PARSE, ...)")
    ap.add_argument("--pocket", type=int, default=None, help="Manually choose fpocket index (1..N); overrides auto picking")
    ap.add_argument("--keep", action="store_true", help="Keep intermediate files")
    args = ap.parse_args()

    receptor = Path(args.receptor).resolve()
    ligand = Path(args.ligand).resolve() if args.ligand else None
    outprefix = Path(args.out).resolve()

    # deps
    fpocket_bin = which_or_die("fpocket", "Install with: conda install -c bioconda -c conda-forge fpocket")
    apbs_bin    = which_or_die("apbs",    "Install with: conda install -c conda-forge apbs")
    pdb2pqr_bin = shutil.which("pdb2pqr") or shutil.which("pdb2pqr30") or shutil.which("pdb2pqr.py")
    if not pdb2pqr_bin:
        raise SystemExit("pdb2pqr not found. Install with: conda install -c conda-forge pdb2pqr")
    print(f"APBS   : {apbs_bin}")
    print(f"PDB2PQR: {pdb2pqr_bin}")

    if not receptor.exists():
        raise SystemExit(f"Receptor not found: {receptor}")
    if ligand and not ligand.exists():
        print(f"WARNING: Ligand file not found; continuing without: {ligand}", file=sys.stderr)
        ligand = None

    workdir = receptor.parent
    print(f"Workdir: {workdir}")

    # 1) fpocket
    print("\n[1/3] Running fpocket ...")
    fp_out_dir = workdir / f"{receptor.stem}_out"
    if fp_out_dir.exists():
        shutil.rmtree(fp_out_dir)
    run([fpocket_bin, "-f", str(receptor)], cwd=workdir)
    if not fp_out_dir.exists():
        raise SystemExit("fpocket output folder not found, something went wrong.")

    # choose pocket
    chosen_idx = args.pocket
    chosen_vert = None
    chosen_vol = None

    pockets_dir = fp_out_dir / "pockets"
    candidates = sorted(pockets_dir.glob("pocket*_vert.pqr"))
    if not candidates:
        print("WARNING: fpocket found no pockets.")

    if chosen_idx is None and candidates:
        ligand_ctr = load_sdf_centroid(str(ligand)) if ligand else None

        # map index -> volume (robust parse)
        vol_map = {}
        info_file = fp_out_dir / "pockets_info.txt"
        if info_file.exists():
            txt = info_file.read_text()
            for m in re.finditer(r"Pocket\s+(\d+).*?Volume\s*:\s*([\d.]+)", txt, re.S | re.I):
                vol_map[int(m.group(1))] = float(m.group(2))

        best_by_vol = (None, -1.0, None)   # (idx, vol, path)
        best_by_vert= (None, -1,   None)   # (idx, nverts, path)
        best = None; best_dist = 1e9

        for p in candidates:
            m = re.search(r"pocket(\d+)_vert\.pqr", p.name)
            if not m: continue
            idx = int(m.group(1))
            vol = vol_map.get(idx, -1.0)
            if vol > best_by_vol[1]:
                best_by_vol = (idx, vol, p)
            # count vertices
            nverts = 0
            try:
                with open(p) as fh:
                    for line in fh:
                        if line.startswith(("ATOM","HETATM")):
                            nverts += 1
            except:
                pass
            if nverts > best_by_vert[1]:
                best_by_vert = (idx, nverts, p)

            if ligand_ctr is not None:
                arr = parse_pqr_xyz(p)
                if arr is not None and len(arr):
                    ctr = arr.mean(axis=0)
                    dist = float(np.linalg.norm(ctr - ligand_ctr))
                    if dist < best_dist:
                        best_dist = dist
                        best = (idx, vol, p)

        if ligand_ctr is not None and best:
            chosen_idx, chosen_vol, chosen_vert = best
            print(f"  auto-picked by proximity: pocket {chosen_idx} (dist ~ {best_dist:.2f} Å, vol ~ {chosen_vol if isinstance(chosen_vol,float) else 'n/a'} Å^3)")
        elif best_by_vol[0] is not None and best_by_vol[1] >= 0:
            chosen_idx, chosen_vol, chosen_vert = best_by_vol
            print(f"  fallback by volume: pocket {chosen_idx} (vol ~ {chosen_vol:.1f} Å^3)")
        elif best_by_vert[0] is not None:
            chosen_idx, _, chosen_vert = best_by_vert
            chosen_vol = vol_map.get(chosen_idx, None)
            print(f"  fallback by vertices: pocket {chosen_idx} (nverts ~ {best_by_vert[1]})")
        else:
            print("WARNING: could not choose a pocket.")
    elif chosen_idx is not None:
        chosen_vert = fp_out_dir / "pockets" / f"pocket{chosen_idx}_vert.pqr"
        info_file = fp_out_dir / "pockets_info.txt"
        if info_file.exists():
            txt = info_file.read_text()
            m = re.search(rf"Pocket\s+{chosen_idx}.*?Volume\s*:\s*([\d.]+)", txt, re.S | re.I)
            if m: chosen_vol = float(m.group(1))
        print(f"  user-picked pocket: {chosen_idx} (vol ~ {chosen_vol if chosen_vol else 'n/a'} Å^3)")

    if not chosen_vert or not chosen_vert.exists():
        print("WARNING: chosen pocket vertices not found; continuing without pocket surface.")

    # 2) PDB2PQR + APBS
    print("\n[2/3] Running PDB2PQR + APBS ...")
    pqr_path = outprefix.with_suffix(".pqr")

    tried = False
    for flags in (
        [f"--ff={args.ff}", f"--with-ph={args.pH}"],
        [f"--ff={args.ff}", "--with-ph", str(args.pH)],
        ["--ff", str(args.ff), "--with-ph", str(args.pH)],
    ):
        try:
            run([pdb2pqr_bin, *flags, str(receptor), str(pqr_path)])
            tried = True
            break
        except SystemExit:
            continue
    if not tried or not pqr_path.exists():
        raise SystemExit("pdb2pqr failed; try a different --ff or check your pdb2pqr version (`--help`).")

    bb = pdb_bbox(str(receptor))
    if not bb:
        raise SystemExit("Failed to parse receptor coordinates.")
    mn, mx = bb
    span = mx - mn
    pad = 20.0
    fg = span + 2*pad          # fine-grid box size (Å)
    cg = fg * 1.4              # coarse-grid box size (Å)
    dime = [161, 161, 161]     # grid points per dimension

    apbs_in = outprefix.with_suffix(".apbs.in")
    dx_out = outprefix.with_name(outprefix.name + "_pot.dx")

    def write_apbs_in(mode="auto"):
        with open(apbs_in, "w") as fh:
            if mode == "auto":
                fh.write(f"""read
  mol pqr {pqr_path.name}
end
elec
  mg-auto
  dime {dime[0]} {dime[1]} {dime[2]}
  cglen {cg[0]:.1f} {cg[1]:.1f} {cg[2]:.1f}
  fglen {fg[0]:.1f} {fg[1]:.1f} {fg[2]:.1f}
  cgcent mol 1
  fgcent mol 1
  mol 1
  lpbe
  bcfl sdh
  pdie 2.0
  sdie 78.5
  srfm smol
  chgm spl2
  sdens 10.0
  srad 1.4
  temp 298.15
  ion charge 1  conc 0.150 radius 2.0
  ion charge -1 conc 0.150 radius 2.0
  write pot dx {dx_out.name.replace('.dx','')}
end
quit
""")
            else:
                # Older APBS fallback (manual grid)
                fh.write(f"""read
  mol pqr {pqr_path.name}
end
elec
  mg-manual
  dime 129 129 129
  grid 0.6 0.6 0.6
  gcent mol 1
  mol 1
  lpbe
  bcfl sdh
  pdie 2.0
  sdie 78.5
  srfm smol
  chgm spl2
  sdens 10.0
  srad 1.4
  temp 298.15
  ion charge 1  conc 0.150 radius 2.0
  ion charge -1 conc 0.150 radius 2.0
  write pot dx {dx_out.name.replace('.dx','')}
end
quit
""")

    # Try mg-auto, then mg-manual fallback if needed
    write_apbs_in(mode="auto")
    try:
        run([apbs_bin, apbs_in.name], cwd=str(outprefix.parent))
    except SystemExit:
        print("\nAPBS mg-auto failed. Retrying with mg-manual...\n")
        write_apbs_in(mode="manual")
        run([apbs_bin, apbs_in.name], cwd=str(outprefix.parent))

    if not dx_out.exists():
        raise SystemExit("APBS finished without producing a .dx map. Check APBS stdout above.")

    # 3) PyMOL: ESP surface + pocket; save PSE/PNG
    print("\n[3/3] Making PyMOL session + PNG ...")
    pse_file = outprefix.with_name(outprefix.name + "_esp.pse")
    png_file = outprefix.with_name(outprefix.name + "_esp.png")
    pymol_session_png(
        str(receptor),
        str(ligand) if ligand else None,
        str(dx_out),
        str(chosen_vert) if chosen_vert and chosen_vert.exists() else None,
        str(pse_file),
        str(png_file)
    )

    # summary
    summary = {
        "receptor": str(receptor),
        "ligand": str(ligand) if ligand else None,
        "fpocket_dir": str(fp_out_dir),
        "chosen_pocket_index": int(chosen_idx) if chosen_idx else None,
        "chosen_pocket_volume_A3": float(chosen_vol) if chosen_vol is not None else None,
        "apbs_dx": str(dx_out),
        "pymol_session": str(pse_file),
        "png": str(png_file),
    }
    with open(outprefix.with_name(outprefix.name + "_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print("\nDone.")
    print(json.dumps(summary, indent=2))

    if not args.keep:
        # keep fpocket outputs and APBS inputs for reproducibility
        pass

if __name__ == "__main__":
    main()
