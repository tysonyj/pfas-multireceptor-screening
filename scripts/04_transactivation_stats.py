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

print(f"\n{'Compound':<10}{'PXR':>9}{'PPAR-a':>10}{'ratio':>8}{'preferred':>11}"
      f"{'p (biol.)':>12}{'p (tech.)':>12}")
print("-" * 72)
print("  p (biol.) uses the three biological replicate means and is the value")
print("  reported in the manuscript. p (tech.) treats all nine measurements as")
print("  independent, which they are not, and overstates significance.")
print("-" * 72)
for c in ['GenX', 'PFECHS', 'F-53B', 'HFPO-TA']:
    sub = d[(d.Compound == c) & (d.Concentration_uM == 100)]
    a = sub[sub.Receptor == 'PXR']
    b = sub[sub.Receptor == 'PPAR-a']

    # technical: all nine measurements
    _, p_tech = stats.ttest_ind(a.Fold_induction, b.Fold_induction, equal_var=False)

    # biological: mean of each biological replicate (replicates 1-3, 4-6, 7-9)
    def bio_means(x):
        x = x.sort_values('Replicate')
        return [x.Fold_induction.iloc[i:i+3].mean() for i in (0, 3, 6)]
    _, p_bio = stats.ttest_ind(bio_means(a), bio_means(b), equal_var=False)

    pref = 'PXR' if a.Fold_induction.mean() > b.Fold_induction.mean() else 'PPAR-a'
    print(f"{c:<10}{a.Fold_induction.mean():>9.2f}{b.Fold_induction.mean():>10.2f}"
          f"{a.Fold_induction.mean()/b.Fold_induction.mean():>8.2f}{pref:>11}"
          f"{p_bio:>12.2e}{p_tech:>12.2e}")
