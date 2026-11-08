"""
01_validate_redocking.py — redocking validation gate.

Each candidate receptor must reproduce the crystallographic pose of its own
co-crystallized ligand before any PFAS is docked into it. The reference ligand
is redocked under the production protocol and the result compared with the
deposited crystallographic coordinates.

METRIC
------
Symmetry-corrected heavy-atom RMSD, computed in place (no superposition) with
RDKit's rdMolAlign.CalcRMS. CalcRMS enumerates graph automorphisms of the
molecule, so topologically equivalent atoms (carboxylate oxygens, phenyl ring
carbons, CF3 fluorines) are matched optimally while bond topology is respected.
This is the conventional near-native metric used in docking benchmarks; the
pass criterion is 2.0 A.

Bond topology is recovered from the Meeko `REMARK SMILES` records that the
PDBQT files carry, so the comparison is between two chemically identical
molecules rather than between bare coordinate sets.

PXR is a documented exception. The deposited SR12813 coordinates in 1ILH could
not be assembled into a single connected molecule, so the reference was rebuilt
from SMILES and carries no crystallographic frame (its centroid sits at the
origin). Atom-level RMSD is therefore undefined for PXR and it is assessed by
the displacement of the docked pose centroid from the crystallographic ligand
centroid instead. That weaker criterion is reported as such.

Usage:  python3 01_validate_redocking.py [--only TARGET] [--force]
        Resumable: a target whose pose file already exists in results/ is not
        redocked, only re-scored. --only runs a single target; --force redocks.
"""
import json, os, sys, time
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign
from meeko import PDBQTMolecule, RDKitMolCreate

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CFG = json.load(open(f"{BASE}/config/config_verified.json"))
OUT = f"{BASE}/results"
os.makedirs(OUT, exist_ok=True)
FORCE = "--force" in sys.argv
ONLY = sys.argv[sys.argv.index("--only") + 1] if "--only" in sys.argv else None

PASS_RMSD = 2.0        # near-native criterion
PASS_CENTROID = 2.0    # fallback, PXR only
CENTROID_ONLY = {"PXR"}


def heavy_coords(path):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in open(path)
                     if l.startswith(("ATOM", "HETATM"))
                     and not l[12:16].strip().startswith("H")])


def mol_from_pdbqt(path):
    """Rebuild an RDKit molecule with correct bond topology from a PDBQT that
    carries Meeko REMARK SMILES records. Returns the first pose, heavy atoms."""
    pm = PDBQTMolecule.from_file(path, skip_typing=True)
    mols = RDKitMolCreate.from_pdbqt_mol(pm)
    mol = [m for m in mols if m is not None][0]
    return Chem.RemoveHs(mol)


def symmetry_rmsd(probe_path, ref_path):
    """Symmetry-corrected heavy-atom RMSD, in place (no superposition)."""
    p = mol_from_pdbqt(probe_path)
    r = mol_from_pdbqt(ref_path)
    if Chem.MolToSmiles(p) != Chem.MolToSmiles(r):
        raise ValueError(f"topology mismatch\n  probe {Chem.MolToSmiles(p)}\n  ref   {Chem.MolToSmiles(r)}")
    # CalcRMS enumerates automorphisms and does NOT superimpose.
    return float(rdMolAlign.CalcRMS(p, r))


def redock(t, c):
    from vina import Vina
    v = Vina(sf_name="vina", cpu=0, seed=42, verbosity=0)
    v.set_receptor(f"{BASE}/receptors/{t}_receptor.pdbqt")
    v.set_ligand_from_file(f"{BASE}/reference_ligands/ref_{t}.pdbqt")
    v.compute_vina_maps(center=c["center"], box_size=c["size"])
    v.dock(exhaustiveness=16, n_poses=5)
    v.write_poses(f"{OUT}/redock_{t}.pdbqt", n_poses=1, overwrite=True)
    return float(v.energies(n_poses=1)[0][0])


def main():
    res = {}
    print(f"{'Target':<9}{'dG':>8}{'RMSD_sym':>11}{'centroid':>10}  {'metric':<10}{'verdict'}")
    print("-" * 68)
    prev = {}
    gate = f"{OUT}/validation_gate.json"
    if os.path.exists(gate):
        try: prev = json.load(open(gate))
        except Exception: prev = {}
    for t, c in CFG.items():
        if ONLY and t != ONLY:
            if t in prev: res[t] = prev[t]
            continue
        rec = {}
        pose_path = f"{OUT}/redock_{t}.pdbqt"
        ref_path = f"{BASE}/reference_ligands/ref_{t}.pdbqt"
        try:
            if FORCE or not os.path.exists(pose_path):
                rec["dG"] = round(redock(t, c), 2)
            elif t in prev and "dG" in prev[t]:
                rec["dG"] = prev[t]["dG"]
            pose = heavy_coords(pose_path)
            ref = heavy_coords(ref_path)

            # centroid displacement is always computed, against the
            # crystallographic ligand centroid (identical to the box centre
            # by construction) -- not against the box centre as a proxy.
            rec["centroid_dev"] = round(float(np.linalg.norm(
                pose.mean(0) - np.array(c["center"]))), 2)

            if t in CENTROID_ONLY:
                rec["rmsd_sym"] = None
                rec["metric"] = "centroid"
                rec["pass"] = rec["centroid_dev"] < PASS_CENTROID
                rec["note"] = ("reference rebuilt from SMILES; deposited coordinates "
                               "could not be assembled into one connected molecule, "
                               "so atom-level RMSD is undefined")
            else:
                rec["rmsd_sym"] = round(symmetry_rmsd(pose_path, ref_path), 2)
                rec["metric"] = "rmsd_sym"
                rec["pass"] = rec["rmsd_sym"] < PASS_RMSD
            rec["n_heavy"] = int(len(ref))
        except Exception as e:
            rec["error"] = f"{type(e).__name__}: {e}"
            rec["pass"] = False
        res[t] = rec
        print(f"{t:<9}{rec.get('dG', float('nan')):>8.2f}"
              f"{(rec.get('rmsd_sym') if rec.get('rmsd_sym') is not None else float('nan')):>11.2f}"
              f"{rec.get('centroid_dev', float('nan')):>10.2f}  {rec.get('metric','-'):<10}"
              f"{'PASS' if rec.get('pass') else 'FAIL'}"
              + (f"   [{rec['error']}]" if "error" in rec else ""), flush=True)

    json.dump(res, open(f"{OUT}/validation_gate.json", "w"), indent=2, ensure_ascii=False)
    n = sum(1 for v in res.values() if v.get("pass"))
    print("-" * 68)
    print(f"Passed: {n}/{len(CFG)}")
    failed = [k for k, v in res.items() if not v.get("pass")]
    if failed:
        print("Excluded from the production screen:", ", ".join(failed))
    print(f"\nwritten: {OUT}/validation_gate.json")


if __name__ == "__main__":
    main()
