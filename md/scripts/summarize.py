"""Summarize the MD metrics for all four systems and compare binding with function."""
import os, glob, json
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
A    = f'{BASE}/analysis'
SYS  = ['PXR_PFECHS','PPARa_PFECHS','PXR_F-53B','PPARa_F-53B']

def read_xvg(p):
    if not os.path.exists(p): return None
    v=[]
    for l in open(p):
        if l.startswith(('#','@')): continue
        try: v.append([float(x) for x in l.split()])
        except: pass
    return np.array(v) if v else None

def plateau(a, col=1, frac=0.3):
    if a is None or len(a)==0: return None
    n=max(1,int(len(a)*frac))
    return float(a[-n:,col].mean())

res={}
for s in SYS:
    d=f'{A}/{s}'
    r={'rmsd_ligand_nm':plateau(read_xvg(f'{d}/rmsd_ligand.xvg')),
       'rmsd_protein_nm':plateau(read_xvg(f'{d}/rmsd_protein.xvg')),
       'mindist_nm':plateau(read_xvg(f'{d}/mindist.xvg')),
       'contacts':plateau(read_xvg(f'{d}/contacts.xvg'))}
    ie=read_xvg(f'{d}/interaction_energy.xvg')
    if ie is not None and ie.shape[1]>=3:
        coul=ie[:,1].mean(); lj=ie[:,2].mean()
        r['coul_kJmol']=float(coul); r['lj_kJmol']=float(lj)
        r['IE_total_kJmol']=float(coul+lj)        # total protein-ligand interaction energy
        r['IE_sd']=float((ie[:,1]+ie[:,2]).std())
    res[s]=r

print("="*78)
print("MD summary (interaction energy = Coul-SR + LJ-SR, protein-ligand; mean over final 30%)")
print("="*78)
print(f"{'System':<15}{'Ligand RMSD':>12}{'Min dist':>11}{'Contacts':>9}{'IE (kJ/mol)':>14}")
print("-"*78)
for s in SYS:
    r=res[s]
    f=lambda k,fmt: (fmt%r[k]) if r.get(k) is not None else '  n/a'
    print(f"{s:<15}{f('rmsd_ligand_nm','%9.3f nm')}{f('mindist_nm','%8.3f nm')}"
          f"{f('contacts','%9.0f')}{f('IE_total_kJmol','%14.1f')}")

print()
print("="*78)
print("Key comparison — does binding predict function?")
print("="*78)
ASSAY={'PFECHS':('PXR',5.29),'F-53B':('PPAR-a',2.88)}
for lg,key in [('PFECHS',('PXR_PFECHS','PPARa_PFECHS')),('F-53B',('PXR_F-53B','PPARa_F-53B'))]:
    a,b=key
    la=res[a].get('IE_total_kJmol'); lb=res[b].get('IE_total_kJmol')
    pref,fold=ASSAY[lg]
    if la is None or lb is None:
        print(f"{lg}: interaction energy not available"); continue
    lie_pref = 'PXR' if la<lb else 'PPAR-a'
    mark = 'match' if lie_pref==pref else 'MISMATCH'
    print(f"{lg:<9} IE: PXR {la:8.1f} / PPAR-a {lb:8.1f}  -> {lie_pref:<7}"
          f" | assay {pref} ({fold:.2f}x)  {mark}")

json.dump(res, open(f'{A}/md_summary.json','w'), indent=2)
print(f"\nSaved: {A}/md_summary.json")
