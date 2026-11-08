"""
04_transactivation_stats.py — statistics for the measured reporter assays.

Replicate structure
-------------------
Each Compound x Receptor x Concentration condition has nine measurements:
three independent biological replicates, each measured in technical triplicate.
Replicates 1-3 belong to biological replicate 1, 4-6 to biological replicate 2,
7-9 to biological replicate 3.

Technical triplicates are averaged within each biological replicate BEFORE
curve fitting and hypothesis testing. The unit of analysis is therefore n = 3
biological replicates, not N = 9 wells.

Scope
-----
The PXR assay (HepG2, native CYP3A4-XREM promoter) and the PPAR-a assay
(COS-1, GAL4-LBD chimera, synthetic UAS reporter) differ in cell background,
receptor construct, promoter architecture and expression level. Response
magnitudes are NOT comparable between them, and a ratio taken within a compound
does not cancel those differences. This script therefore performs only
within-receptor analysis and makes no receptor-preference assignment.

Usage:  python3 04_transactivation_stats.py
"""
import os
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
d = pd.read_csv(f'{BASE}/data/experimental/transactivation_source_data.csv')

# replicate 1-9 -> biological replicate 1-3
d['Bio'] = ((d.Replicate - 1) // 3) + 1
d['Tech'] = ((d.Replicate - 1) % 3) + 1

COMPOUNDS = ['GenX', 'PFECHS', 'F-53B', 'HFPO-TA']
RECEPTORS = ['PXR', 'PPAR-a']


def fpl(x, bottom, top, ec50, hill):
    return bottom + (top - bottom) / (1 + (ec50 / x) ** hill)


def holm(pvals):
    """Holm-Bonferroni adjusted p-values, order preserved."""
    p = np.asarray(pvals, float)
    m = len(p)
    order = np.argsort(p)
    adj = np.empty(m)
    running = 0.0
    for rank, idx in enumerate(order):
        val = (m - rank) * p[idx]
        running = max(running, val)
        adj[idx] = min(running, 1.0)
    return adj


# ----------------------------------------------------------------------
# 1. Replicate variability: technical and biological RSD
# ----------------------------------------------------------------------
print("=" * 78)
print("1. REPLICATE VARIABILITY (relative standard deviation, %)")
print("=" * 78)
print("  Technical RSD : SD / mean x 100 among the three technical wells,")
print("                  computed within each biological replicate and then")
print("                  averaged over the three biological replicates.")
print("  Biological RSD: SD / mean x 100 among the three biological replicate")
print("                  means.")
print("-" * 78)
print(f"{'Compound':<10}{'Receptor':<9}{'Conc (uM)':>10}{'mean FI':>10}"
      f"{'tech RSD':>10}{'biol RSD':>10}")
print("-" * 78)

rsd_rows = []
for c in COMPOUNDS:
    for r in RECEPTORS:
        sub = d[(d.Compound == c) & (d.Receptor == r)]
        for conc, g in sub.groupby('Concentration_uM'):
            bio_means, tech_rsds = [], []
            for b, gb in g.groupby('Bio'):
                v = gb.Fold_induction.values
                bio_means.append(v.mean())
                tech_rsds.append(100 * v.std(ddof=1) / v.mean())
            bio_means = np.array(bio_means)
            tech_rsd = float(np.mean(tech_rsds))
            biol_rsd = 100 * bio_means.std(ddof=1) / bio_means.mean()
            rsd_rows.append(dict(Compound=c, Receptor=r, Conc=conc,
                                 mean=bio_means.mean(),
                                 tech_rsd=tech_rsd, biol_rsd=biol_rsd))
            print(f"{c:<10}{r:<9}{conc:>10.1f}{bio_means.mean():>10.2f}"
                  f"{tech_rsd:>10.1f}{biol_rsd:>10.1f}")

rsd = pd.DataFrame(rsd_rows)
print("-" * 78)
print(f"{'ALL CONDITIONS':<29}{'':>10}{'':>10}"
      f"{rsd.tech_rsd.median():>10.1f}{rsd.biol_rsd.median():>10.1f}   (median)")
print(f"{'':<29}{'':>10}{'':>10}"
      f"{rsd.tech_rsd.max():>10.1f}{rsd.biol_rsd.max():>10.1f}   (max)")

rsd.to_csv(f'{BASE}/data/experimental/replicate_rsd.csv', index=False)
print(f"\nwritten: data/experimental/replicate_rsd.csv")


# ----------------------------------------------------------------------
# 2. Four-parameter logistic fits (within receptor)
# ----------------------------------------------------------------------
print("\n" + "=" * 78)
print("2. FOUR-PARAMETER LOGISTIC FITS (fitted to biological replicate means)")
print("=" * 78)
print(f"{'Compound':<10}{'Receptor':<9}{'EC50 (uM)':>22}{'Emax':>20}{'Hill':>16}")
print("-" * 78)

fit_rows = []
for c in COMPOUNDS:
    for r in RECEPTORS:
        sub = d[(d.Compound == c) & (d.Receptor == r)]
        # technical -> biological -> condition mean
        bio = sub.groupby(['Concentration_uM', 'Bio']).Fold_induction.mean()
        g = bio.groupby('Concentration_uM').mean().sort_index()
        xv, yv = g.index.values.astype(float), g.values
        try:
            p, cov = curve_fit(fpl, xv, yv,
                               p0=[1, yv.max() * 1.2, 20, 1.5], maxfev=200000)
            se = np.sqrt(np.diag(cov))
            dof = max(len(xv) - 4, 1)
            t = stats.t.ppf(0.975, dof)
            ci = t * se
            fit_rows.append(dict(Compound=c, Receptor=r,
                                 ec50=p[2], ec50_se=se[2], ec50_lo=p[2] - ci[2],
                                 ec50_hi=p[2] + ci[2],
                                 emax=p[1], emax_se=se[1], emax_lo=p[1] - ci[1],
                                 emax_hi=p[1] + ci[1],
                                 hill=p[3], hill_se=se[3]))
            print(f"{c:<10}{r:<9}"
                  f"{p[2]:>8.1f} [{p[2]-ci[2]:>6.1f},{p[2]+ci[2]:>6.1f}]"
                  f"{p[1]:>8.1f} [{p[1]-ci[1]:>4.1f},{p[1]+ci[1]:>4.1f}]"
                  f"{p[3]:>8.2f} +-{se[3]:<5.2f}")
        except Exception as e:
            print(f"{c:<10}{r:<9}{'fit failed':>22}   ({e})")

pd.DataFrame(fit_rows).to_csv(
    f'{BASE}/data/experimental/dose_response_fits.csv', index=False)
print("\nBrackets are 95% confidence intervals from the fit covariance matrix")
print("(t distribution, 8 concentrations - 4 parameters = 4 dof).")
print("written: data/experimental/dose_response_fits.csv")


# ----------------------------------------------------------------------
# 3. Within-receptor activation tests, Holm-corrected
# ----------------------------------------------------------------------
print("\n" + "=" * 78)
print("3. WITHIN-RECEPTOR ACTIVATION (100 uM vs lowest tested concentration)")
print("=" * 78)
print("  Paired t-test on the three biological replicate means. Eight tests")
print("  (4 compounds x 2 receptors) corrected by Holm-Bonferroni.")
print("  No comparison is made between the two receptors.")
print("-" * 78)

labels, raw_p, effect = [], [], []
for c in COMPOUNDS:
    for r in RECEPTORS:
        sub = d[(d.Compound == c) & (d.Receptor == r)]
        lo_conc = sub.Concentration_uM.min()
        hi = sub[sub.Concentration_uM == 100].groupby('Bio').Fold_induction.mean()
        lo = sub[sub.Concentration_uM == lo_conc].groupby('Bio').Fold_induction.mean()
        t_stat, p = stats.ttest_rel(hi.values, lo.values)
        labels.append((c, r))
        raw_p.append(p)
        effect.append(hi.mean() / lo.mean())

adj_p = holm(raw_p)
bonf_p = np.minimum(np.asarray(raw_p) * len(raw_p), 1.0)

print(f"{'Compound':<10}{'Receptor':<9}{'FI(100)/FI(low)':>17}"
      f"{'p raw':>12}{'p Holm':>12}{'p Bonf':>12}")
print("-" * 78)
for (c, r), e, p0, ph, pb in zip(labels, effect, raw_p, adj_p, bonf_p):
    star = '' if ph < 0.05 else '   n.s.'
    print(f"{c:<10}{r:<9}{e:>17.2f}{p0:>12.2e}{ph:>12.2e}{pb:>12.2e}{star}")

print("-" * 78)
print(f"Surviving Holm correction at alpha = 0.05: "
      f"{int((adj_p < 0.05).sum())} of {len(adj_p)}")
print(f"Surviving Bonferroni correction at alpha = 0.05: "
      f"{int((bonf_p < 0.05).sum())} of {len(bonf_p)}")


# ----------------------------------------------------------------------
# 4. LaTeX table for the Supporting Information
# ----------------------------------------------------------------------
out = [r'\begin{tabular}{llrrrr}', r'\toprule',
       r'Compound & Receptor & Conc.\ ($\mu$M) & Mean FI & '
       r'Technical RSD (\%) & Biological RSD (\%) \\', r'\midrule']
for c in COMPOUNDS:
    for r in RECEPTORS:
        s = rsd[(rsd.Compound == c) & (rsd.Receptor == r)].sort_values('Conc')
        rec = r.replace('PPAR-a', r'PPAR-\textalpha')
        for i, row in enumerate(s.itertuples()):
            c_cell = c if i == 0 else ''
            r_cell = rec if i == 0 else ''
            out.append(f'{c_cell} & {r_cell} & {row.Conc:.1f} & {row.mean:.2f} & '
                       f'{row.tech_rsd:.1f} & {row.biol_rsd:.1f} \\\\')
        out.append(r'\midrule')
out[-1] = r'\bottomrule'
out.append(r'\end{tabular}')

with open(f'{BASE}/data/experimental/tableS7_rsd.tex', 'w') as f:
    f.write('\n'.join(out) + '\n')
print("\nwritten: data/experimental/tableS7_rsd.tex  (LaTeX table for the SI)")
