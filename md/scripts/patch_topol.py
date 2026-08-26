"""Add the ligand itp #include and a [ molecules ] entry to topol.top."""
import argparse, os, re
p=argparse.ArgumentParser()
p.add_argument('--topol',required=True); p.add_argument('--lig-itp',required=True)
p.add_argument('--lig-name',default='LIG')
a=p.parse_args()
top=open(a.topol).read()
if a.lig_name in top and '#include' in top and os.path.basename(a.lig_itp) in top:
    print("  topol.top already patched"); raise SystemExit

itp=open(a.lig_itp).read()
# Split out atomtypes (must come immediately after the forcefield include)
m=re.search(r'\[\s*atomtypes\s*\][^\[]*', itp)
atomtypes = m.group(0) if m else ''
rest = itp.replace(atomtypes,'') if atomtypes else itp
d=os.path.dirname(a.topol) or '.'
open(f'{d}/lig_atomtypes.itp','w').write(atomtypes)
open(f'{d}/lig.itp','w').write(rest)

lines=top.splitlines(); out=[]; done_at=False
for l in lines:
    out.append(l)
    if not done_at and l.strip().startswith('#include') and 'forcefield.itp' in l:
        out += ['', '; ligand atomtypes', '#include "lig_atomtypes.itp"',
                '', '; ligand topology', '#include "lig.itp"', '']
        done_at=True
top="\n".join(out)
top=re.sub(r'(\[\s*molecules\s*\][^\n]*\n(?:[^\n]*\n)*?)(\s*$)',
           lambda m: m.group(1)+f'{a.lig_name}                 1\n', top, count=1)
if f'\n{a.lig_name}' not in top.split('[ molecules ]')[-1]:
    top = top.rstrip()+f'\n{a.lig_name}                 1\n'
open(a.topol,'w').write(top)
print(f"  topol.top patched ({a.lig_name} added)")
