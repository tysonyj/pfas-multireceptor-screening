"""
03_statistics.py — correlation analysis and leverage diagnostics for the
docking results.
Usage:  python3 03_statistics.py
"""
import os
import numpy as np, pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
d = pd.read_csv(f'{BASE}/data/docking/analysis_input.csv')
x, y, z = d.s_bind.values, d.hazard_H.values, d.nF.values.astype(float)

def partial_spearman(x, y, z):
    rxy = stats.spearmanr(x, y)[0]
    rxz = stats.spearmanr(x, z)[0]
    ryz = stats.spearmanr(y, z)[0]
    return (rxy - rxz*ryz) / np.sqrt((1 - rxz**2) * (1 - ryz**2)), rxy, rxz, ryz

pr, rxy, rxz, ryz = partial_spearman(x, y, z)
n  = len(d); dfree = n - 3
pv = 2 * (1 - stats.t.cdf(abs(pr * np.sqrt(dfree / (1 - pr**2))), dfree))
_, pxy = stats.spearmanr(x, y)

rng, bs = np.random.default_rng(42), []
for _ in range(5000):
    i = rng.integers(0, n, n)
    try:
        v = partial_spearman(x[i], y[i], z[i])[0]
        if np.isfinite(v): bs.append(v)
    except Exception: pass
lo, hi = np.percentile(bs, [2.5, 97.5])

K  = 10
ov = len(set(d.nlargest(K, 's_bind').pfas) & set(d.nlargest(K, 'hazard_H').pfas))
orr, pf = stats.fisher_exact([[ov, K - ov], [K - ov, n - 2*K + ov]])

print(f"n = {n} compounds\n")
print(f"Spearman rho (s_bind vs hazard_H)   = {rxy:+.3f}   p = {pxy:.2e}")
print(f"  r(s_bind, nF)                     = {rxz:+.3f}")
print(f"  r(hazard_H, nF)                   = {ryz:+.3f}")
print(f"Partial rho controlling for nF      = {pr:+.3f}   p = {pv:.4f}")
print(f"  bootstrap 95% CI                  = [{lo:+.3f}, {hi:+.3f}]")
print(f"Fisher top-{K} overlap               = {ov}/{K}   OR = {orr:.2f}   p = {pf:.3f}")

print("\nLeave-one-out leverage on partial rho (top 6):")
rows = []
for i in range(n):
    sub = d.drop(d.index[i])
    v = partial_spearman(sub.s_bind.values, sub.hazard_H.values,
                         sub.nF.values.astype(float))[0]
    rows.append((d.pfas.iloc[i], v, abs(v - pr)))
for nm, v, dl in sorted(rows, key=lambda r: -r[2])[:6]:
    print(f"  without {nm:<10} partial rho = {v:+.3f}   (change {dl:.3f})")
