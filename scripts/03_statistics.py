"""
03_statistics.py — association between cross-target binding strength and
regulatory priority, controlling for perfluorocarbon chain length.

The dependent variable is `reg_score`, a count of binding regulatory
instruments in which each compound is individually named or falls within an
explicit homologue range (see data/hazard/SOURCES.md for the scoring rule,
the instruments, and the compilation date). It is documentary rather than
estimated, and it is not derived from molecular structure -- which is why it
can be correlated against docking scores without circularity.

Usage:  python3 03_statistics.py
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
d = pd.read_csv(f"{BASE}/data/docking/analysis_input.csv")
n = len(d)


def partial_spearman(x, y, z):
    """Partial Spearman of x,y controlling z, via the three-variable formula."""
    rxy, _ = stats.spearmanr(x, y)
    rxz, _ = stats.spearmanr(x, z)
    ryz, _ = stats.spearmanr(y, z)
    return (rxy - rxz * ryz) / np.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))


print(f"n = {n} compounds (all docked compounds retained; reg_score is defined for every one)\n")
print("reg_score distribution:",
      ", ".join(f"{k}:{v}" for k, v in sorted(d.reg_score.value_counts().items())))

r_sh, p_sh = stats.spearmanr(d.s_bind, d.reg_score)
r_sn, p_sn = stats.spearmanr(d.s_bind, d.nF)
r_hn, p_hn = stats.spearmanr(d.reg_score, d.nF)
pr = partial_spearman(d.s_bind, d.reg_score, d.nF)
dof = n - 3
t = pr * np.sqrt(dof / (1 - pr ** 2))
p_pr = 2 * (1 - stats.t.cdf(abs(t), dof))

print("\n--- rank correlations ---")
print(f"  s_bind   vs reg_score   rho = {r_sh:+.3f}   p = {p_sh:.3g}")
print(f"  s_bind   vs nF          rho = {r_sn:+.3f}   p = {p_sn:.3g}")
print(f"  reg_score vs nF         rho = {r_hn:+.3f}   p = {p_hn:.3g}"
      f"   <- the dependent variable is not a chain-length proxy")
print(f"  partial (control nF)    rho = {pr:+.3f}   p = {p_pr:.3f}")

# bootstrap on the partial correlation
rng = np.random.default_rng(42)
bs = []
for _ in range(5000):
    s = d.iloc[rng.integers(0, n, n)]
    if s.reg_score.nunique() < 2 or s.nF.nunique() < 2:
        continue
    v = partial_spearman(s.s_bind, s.reg_score, s.nF)
    if np.isfinite(v):
        bs.append(v)
lo, hi = np.percentile(bs, [2.5, 97.5])
print(f"  bootstrap 95% CI        [{lo:+.3f}, {hi:+.3f}]   "
      f"({100*np.mean(np.array(bs) < 0):.0f}% of resamples below zero, {len(bs)} valid)")

# leave-one-out leverage
p0 = partial_spearman(d.s_bind, d.reg_score, d.nF)
loo = sorted((partial_spearman(*d.drop(i)[["s_bind", "reg_score", "nF"]].T.values), d.loc[i, "pfas"])
             for i in d.index)
print(f"\n--- leave-one-out leverage on the partial correlation ---")
print(f"  full panel {p0:+.3f};  range {loo[0][0]:+.3f} ({loo[0][1]}) to {loo[-1][0]:+.3f} ({loo[-1][1]})")

# listed vs unlisted
reg, non = d[d.reg_score > 0], d[d.reg_score == 0]
U, p_u = stats.mannwhitneyu(reg.s_bind, non.s_bind, alternative="two-sided")
print(f"\n--- compounds under binding regulation vs those under none ---")
print(f"  listed    n = {len(reg):2d}   s_bind median {reg.s_bind.median():.2f}   nF mean {reg.nF.mean():.1f}")
print(f"  unlisted  n = {len(non):2d}   s_bind median {non.s_bind.median():.2f}   nF mean {non.nF.mean():.1f}")
print(f"  Mann-Whitney U = {U:.0f}, p = {p_u:.3f}")

top10 = d.nlargest(10, "s_bind")
print(f"\n--- the ten strongest binders and their regulatory status ---")
for _, r in top10.iterrows():
    print(f"  {r.pfas:10s} s_bind {r.s_bind:5.2f}   reg_score {int(r.reg_score)}")
print(f"  {int((top10.reg_score == 0).sum())} of the ten strongest binders are named in no binding instrument.")
