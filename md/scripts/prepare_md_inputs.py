"""
prepare_md_inputs.py — regenerate the MD inputs from the anionic docking poses

Usage:
    python3 prepare_md_inputs.py --poses ~/pfas_recovery/results/poses_anionic \
                                 --receptors ~/pfas_recovery/receptors

Replaces the existing systems/*/ligand.mol2 (neutral form) with the anionic
structure.
"""
import argparse, os, shutil, json, sys

ap = argparse.ArgumentParser()
ap.add_argument('--poses',     required=True, help='directory of anionic docking poses')
ap.add_argument('--receptors', required=True, help='directory of receptor pdbqt files')
a = ap.parse_args()

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
POSE = os.path.expanduser(a.poses)
REC  = os.path.expanduser(a.receptors)

try:
    from openbabel import pybel
except ImportError:
    sys.exit("openbabel is required:  conda install -c conda-forge openbabel")

# system name -> (docking target name, ligand name)
SYS = {
    'PXR_PFECHS':   ('PXR',    'PFECHS'),
    'PXR_F-53B':    ('PXR',    'F-53B'),
    'PPARa_PFECHS': ('PPAR-a', 'PFECHS'),
    'PPARa_F-53B':  ('PPAR-a', 'F-53B'),
}

print(f"pose directory     : {POSE}")
print(f"receptor directory : {REC}\n")

ok = 0
for name, (rec, lg) in SYS.items():
    d = f'{BASE}/systems/{name}'
    os.makedirs(d, exist_ok=True)

    src = f'{POSE}/{rec}_{lg}.pdbqt'
    if not os.path.exists(src):
        print(f'{name:<14} pose not found: {src}')
        continue

    # --- Ligand: extract the top-ranked pose only ---
    lines, keep = [], True
    for l in open(src):
        if l.startswith('MODEL') and len(lines) > 0:
            break
        lines.append(l)
    open(f'{d}/_pose1.pdbqt', 'w').writelines(lines)

    m = next(pybel.readfile('pdbqt', f'{d}/_pose1.pdbqt'))
    m.OBMol.AddHydrogens(False, True, 7.4)
    m.write('pdb',  f'{d}/ligand.pdb',  overwrite=True)
    m.write('mol2', f'{d}/ligand.mol2', overwrite=True)
    os.remove(f'{d}/_pose1.pdbqt')

    # --- Receptor: pdbqt -> pdb ---
    rsrc = f'{REC}/{rec}_receptor.pdbqt'
    if os.path.exists(rsrc):
        r = next(pybel.readfile('pdbqt', rsrc))
        r.write('pdb', f'{d}/protein.pdb', overwrite=True)
    else:
        print(f'  Warning: receptor not found {rsrc}')

    dG = None
    for l in open(src):
        if l.startswith('REMARK VINA RESULT'):
            dG = float(l.split()[3]); break

    na = sum(1 for l in open(f'{d}/ligand.pdb') if l.startswith(('ATOM', 'HETATM')))
    json.dump({'system': name, 'receptor': rec, 'ligand': lg,
               'vina_dG': dG, 'ligand_atoms_with_H': na,
               'net_charge': -1,
               'protonation': 'sulfonate deprotonated at pH 7.4 (anionic)'},
              open(f'{d}/info.json', 'w'), indent=2)
    print(f'{name:<14} dG={dG:7.3f}  ligand {na:3d} atoms (incl. H)  -> {d}')
    ok += 1

print(f'\nDone: {ok}/{len(SYS)} systems')
print('Next: bash build_system.sh PXR_PFECHS')
