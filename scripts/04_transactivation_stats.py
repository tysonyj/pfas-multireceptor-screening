"""
04_transactivation_stats.py — statistics for the measured reporter assays.
4PL fits (EC50, Emax, Hill), receptor preference ratios, Welch tests.
Usage:  python3 04_transactivation_stats.py
"""
import os
import numpy as np, pandas as pd
from scipy.optimize import curve_fit
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
d = pd.read_csv(f'{BASE}/data/experimental/transactivation_source_data.csv')

def fpl(x, bottom, top, ec50, hill):
    return bottom + (top - bottom) / (1 + (ec50 / x) ** hill)

print(f"{'Compound':<10}{'Receptor':<9}{'EC50 (uM)':>16}{'Emax':>15}{'Hill':>8}")
print("-" * 60)
fits = {}
for (c, r), sub in d.groupby(['Compound', 'Receptor']):
    g = sub.groupby('Concentration_uM').Fold_induction.mean().sort_index()
    xv, yv = g.index.values.astype(float), g.values
    try:
        p, cov = curve_fit(fpl, xv, yv, p0=[1, yv.max()*1.2, 20, 1.5], maxfev=80000)
        se = np.sqrt(np.diag(cov))
        fits[(c, r)] = (p, se)
        print(f"{c:<10}{r:<9}{p[2]:>9.1f} ± {se[2]:<5.1f}{p[1]:>9.1f} ± {se[1]:<4.1f}{p[3]:>8.2f}")
    except Exception:
        print(f"{c:<10}{r:<9}{'fit failed':>16}")

print(f"\n{'Compound':<10}{'PXR':>9}{'PPAR-a':>10}{'ratio':>9}{'preferred':>12}{'p (Welch)':>13}")
print("-" * 63)
for c in ['GenX', 'PFECHS', 'F-53B', 'HFPO-TA']:
    a = d[(d.Compound==c) & (d.Receptor=='PXR')    & (d.Concentration_uM==100)].Fold_induction
    b = d[(d.Compound==c) & (d.Receptor=='PPAR-a') & (d.Concentration_uM==100)].Fold_induction
    t, p = stats.ttest_ind(a, b, equal_var=False)
    pref = 'PXR' if a.mean() > b.mean() else 'PPAR-a'
    print(f"{c:<10}{a.mean():>9.2f}{b.mean():>10.2f}{a.mean()/b.mean():>9.2f}{pref:>12}{p:>13.2e}")
