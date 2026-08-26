"""Merge a protein .gro and a ligand .gro into a single complex."""
import argparse
p=argparse.ArgumentParser()
p.add_argument('--protein',required=True); p.add_argument('--ligand',required=True)
p.add_argument('--out',required=True)
a=p.parse_args()
prot=open(a.protein).read().splitlines()
lig =open(a.ligand ).read().splitlines()
np_=int(prot[1]); nl=int(lig[1])
body = prot[2:2+np_] + lig[2:2+nl]
box  = prot[2+np_]
with open(a.out,'w') as f:
    f.write("Protein-ligand complex\n")
    f.write(f"{np_+nl:5d}\n")
    f.write("\n".join(body)+"\n")
    f.write(box+"\n")
print(f"  merged: protein {np_} + ligand {nl} = {np_+nl} atoms -> {a.out}")
