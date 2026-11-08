"""
05_make_figures_archtox.py — regenerate the four main figures for the
Archives of Toxicology submission.

Changes relative to the earlier figure set:
  Fig 1  correct PDB codes (previously every row read "1ILH"); grey footnote
         removed (it belongs in the caption); "2.0 Å criterion" label moved
         inside the axes.
  Fig 2  larger fonts; partial-correlation p recomputed with dof = n-3
         (the earlier figure used n-2 and printed p = 0.166).
  Fig 3  larger fonts; ADONA dropped so the panel matches the caption, which
         refers to the four emerging replacement analogs.
  Fig 4  panel C deleted and the cross-assay ratio brackets removed from
         panel B. Both reported a receptor-preference comparison that the
         manuscript withdraws: the two assays use different cell backgrounds
         and reporter architectures and are not comparable on magnitude.
         Fitted parameters now live in Table S7.

Usage:  python3 05_make_figures_archtox.py [outdir]
"""
import os, sys, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
from scipy.optimize import curve_fit

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = sys.argv[1] if len(sys.argv) > 1 else BASE

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 15,
    "axes.labelsize": 16,
    "axes.titlesize": 17,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 14,
    "axes.linewidth": 1.1,
    "savefig.dpi": 400,
    "savefig.bbox": "tight",
})
PANEL = dict(fontsize=22, fontweight="bold", va="top", ha="right")


# ----------------------------------------------------------------------
# Figure 1 — redocking validation
# ----------------------------------------------------------------------
def fig1_validation():
    # Target, PDB, reference ligand, symmetry-corrected RMSD, centroid-only flag
    rows = [
        ("HSA",      "2BXD", "warfarin",      0.45, False),
        ("PPAR-$\\alpha$", "1K7L", "GW409544", 0.52, False),
        ("TR-$\\alpha$",   "2H79", "T3",       0.70, False),
        ("ER-$\\alpha$",   "1GWR", "estradiol", 0.70, False),
        ("TTR",      "1IE4", "thyroxine",     1.19, False),
        ("PPAR-$\\gamma$", "1FM6", "rosiglitazone", 1.97, False),
        ("PXR",      "1ILH", "SR12813",       0.94, True),
        ("L-FABP",   "3B2H", "palmitate",    12.68, False),
    ]
    order = sorted(rows, key=lambda r: (r[3] > 2.0, r[3]))
    labels = [f"{t}\n{lig} ({pdb})" for t, pdb, lig, _, _ in order]
    vals   = [r[3] for r in order]
    passed = [r[3] <= 2.0 for r in order]
    centro = [r[4] for r in order]

    fig, ax = plt.subplots(figsize=(10.5, 7.2))
    y = np.arange(len(order))[::-1]
    for yi, v, ok, c in zip(y, vals, passed, centro):
        ax.barh(yi, v, height=0.62,
                color="#2e7d32" if ok else "#c62828",
                hatch="//" if c else None,
                edgecolor="white" if c else "none", linewidth=0)
        ax.text(v + 0.09, yi, f"{v:.2f}" + ("*" if c else ""),
                va="center", ha="left", fontsize=15)

    ax.axvline(2.0, ls="--", lw=2.0, color="#333333", zorder=3)
    ax.text(2.0 - 0.12, y[0] + 0.55, "2.0 Å criterion", rotation=90,
            va="top", ha="right", fontsize=14, color="#333333")

    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=14)
    ax.set_xlabel("Redocking RMSD (Å)")
    ax.set_xlim(0, 14.2)
    ax.set_ylim(y[-1] - 0.75, y[0] + 0.85)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(os.path.join(OUT, "figAT1_validation.png"))
    plt.close(fig)
    print("figAT1_validation.png")


# ----------------------------------------------------------------------
# Figure 2 — binding strength vs hazard, chain-length confound
# ----------------------------------------------------------------------
CLS = {"PFCA": "#1f77b4", "PFSA": "#d62728", "Ether-linked": "#9467bd",
       "Cyclic": "#2ca02c", "Fluorotelomer": "#ff7f0e", "Sulfonamide": "#7f7f7f"}

def _klass(name):
    if name in ("GenX", "HFPO-TA", "ADONA", "F-53B"): return "Ether-linked"
    if name == "PFECHS": return "Cyclic"
    if "FOSA" in name: return "Sulfonamide"
    if "FT" in name: return "Fluorotelomer"
    if name.endswith("S"): return "PFSA"
    return "PFCA"

def fig2_concordance():
    """Binding strength vs regulatory priority, and the chain-length control."""
    d = pd.read_csv(f"{BASE}/data/docking/analysis_input.csv")
    d["cls"] = d.pfas.map(_klass)
    n = len(d)
    r_sh, p_sh = stats.spearmanr(d.s_bind, d.reg_score)
    r_sn, _ = stats.spearmanr(d.s_bind, d.nF)
    r_hn, p_hn = stats.spearmanr(d.reg_score, d.nF)
    pr = (r_sh - r_sn * r_hn) / np.sqrt((1 - r_sn**2) * (1 - r_hn**2))
    dof = n - 3
    p_pr = 2 * (1 - stats.t.cdf(abs(pr * np.sqrt(dof / (1 - pr**2))), dof))

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.9))
    rng = np.random.default_rng(7)
    jit = lambda k: rng.uniform(-0.13, 0.13, k)

    # --- A: binding strength vs regulatory score -------------------------
    for cls, g in d.groupby("cls"):
        axes[0].scatter(g.s_bind, g.reg_score + jit(len(g)), s=95, color=CLS[cls],
                        label=cls, edgecolor="white", linewidth=0.8, zorder=3)
    axes[0].set_xlabel("Cross-target binding strength, $s_{bind}$")
    axes[0].set_ylabel("Regulatory instruments naming\nthe compound")
    axes[0].set_yticks(range(5))
    axes[0].set_ylim(-0.55, 4.9)
    axes[0].text(0.03, 0.97, f"$\\rho$ = {r_sh:+.3f}   $p$ = {p_sh:.2f}  (n.s.)",
                 transform=axes[0].transAxes, va="top", fontsize=14)

    # --- B: chain-length control, on one common percentile-rank axis -----
    pct = lambda v: 100 * (stats.rankdata(v) - 0.5) / len(v)
    axes[1].scatter(d.nF, pct(d.s_bind), s=95, color="#333333", edgecolor="white",
                    linewidth=0.8, zorder=3,
                    label=f"$s_{{bind}}$   $\\rho$ = {r_sn:+.3f}")
    axes[1].scatter(d.nF, pct(d.reg_score), s=95, color="#c0392b", marker="s",
                    edgecolor="white", linewidth=0.8, zorder=3,
                    label=f"regulatory score   $\\rho$ = {r_hn:+.3f} (n.s.)")
    axes[1].set_xlabel("Number of fluorine atoms, $n_F$")
    axes[1].set_ylabel("Percentile rank within panel")
    axes[1].set_ylim(-6, 118)
    axes[1].legend(frameon=False, loc="upper left", fontsize=13,
                   handletextpad=0.3, borderpad=0.1)

    # --- C: rank residuals after removing nF -----------------------------
    rk = lambda v: stats.rankdata(v)
    rx = rk(d.s_bind) - np.polyval(np.polyfit(rk(d.nF), rk(d.s_bind), 1), rk(d.nF))
    ry = rk(d.reg_score) - np.polyval(np.polyfit(rk(d.nF), rk(d.reg_score), 1), rk(d.nF))
    for cls, g in d.groupby("cls"):
        i = d.index.get_indexer(g.index)
        axes[2].scatter(rx[i], ry[i], s=95, color=CLS[cls], edgecolor="white",
                        linewidth=0.8, zorder=3)
    axes[2].axhline(0, color="#aaaaaa", lw=1); axes[2].axvline(0, color="#aaaaaa", lw=1)
    axes[2].set_xlabel("$s_{bind}$ rank residual  ($n_F$ removed)")
    axes[2].set_ylabel("Regulatory rank residual")
    ylo, yhi = axes[2].get_ylim(); axes[2].set_ylim(ylo, yhi + 0.26 * (yhi - ylo))
    axes[2].text(0.03, 0.97, f"partial $\\rho$ = {pr:+.3f}   $p$ = {p_pr:.2f}  (n.s.)",
                 transform=axes[2].transAxes, va="top", fontsize=14)

    for a, L in zip(axes, "ABC"):
        a.spines[["top", "right"]].set_visible(False)
        a.text(-0.14, 1.11, L, transform=a.transAxes, **PANEL)

    # one shared class legend below the row, so nothing overlaps the data
    handles = [Line2D([], [], marker="o", ls="", ms=11, mec="white", mew=0.8,
                      color=CLS[c], label=c) for c in
               ["PFCA", "PFSA", "Ether-linked", "Cyclic", "Fluorotelomer", "Sulfonamide"]]
    fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False,
               fontsize=14, bbox_to_anchor=(0.5, -0.015), columnspacing=1.6,
               handletextpad=0.35)
    fig.tight_layout(rect=[0, 0.075, 1, 1])
    fig.savefig(os.path.join(OUT, "figAT2_concordance.png"))
    plt.close(fig)
    print(f"figAT2_concordance.png   (rho={r_sh:+.3f} p={p_sh:.3f}; partial={pr:+.3f} p={p_pr:.3f}; rho(reg,nF)={r_hn:+.3f})")


# ----------------------------------------------------------------------
# Figure 3 — target profiles of the four emerging replacement analogs
# ----------------------------------------------------------------------
FOUR = ["GenX", "PFECHS", "F-53B", "HFPO-TA"]
TGT = ["PXR", "PPAR-a", "PPAR-g", "ER-a", "TR-a", "TTR", "HSA"]
TGT_LAB = ["PXR", "PPAR-$\\alpha$", "PPAR-$\\gamma$", "ER-$\\alpha$", "TR-$\\alpha$", "TTR", "HSA"]

def fig3_selectivity():
    dg = pd.read_csv(f"{BASE}/data/docking/dG_matrix.csv", index_col=0).loc[FOUR, TGT]
    fig, axes = plt.subplots(1, 2, figsize=(16, 6.4),
                             gridspec_kw={"width_ratios": [1.15, 1]})
    w = 0.2
    x = np.arange(len(TGT))
    cols = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"]
    for i, (cmp_, c) in enumerate(zip(FOUR, cols)):
        axes[0].bar(x + (i - 1.5) * w, dg.loc[cmp_], w, label=cmp_, color=c)
    axes[0].set_xticks(x); axes[0].set_xticklabels(TGT_LAB, fontsize=14)
    axes[0].set_ylabel("$\\Delta G$ (kcal mol$^{-1}$)")
    axes[0].set_ylim(-11, 0)
    axes[0].legend(frameon=False, ncol=4, loc="lower center", fontsize=14,
                   bbox_to_anchor=(0.5, -0.015), columnspacing=1.2)
    axes[0].spines[["top", "right"]].set_visible(False)

    z = dg.sub(dg.mean(axis=1), axis=0).div(dg.std(axis=1, ddof=0), axis=0)
    im = axes[1].imshow(z.values, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    axes[1].set_xticks(range(len(TGT))); axes[1].set_xticklabels(TGT_LAB, rotation=35, ha="right", fontsize=14)
    axes[1].set_yticks(range(len(FOUR))); axes[1].set_yticklabels(FOUR, fontsize=15)
    for i in range(len(FOUR)):
        for j in range(len(TGT)):
            v = z.values[i, j]
            axes[1].text(j, i, f"{v:+.1f}", ha="center", va="center", fontsize=14,
                         color="white" if abs(v) > 1.3 else "#222222")
    cb = fig.colorbar(im, ax=axes[1], fraction=0.045, pad=0.03)
    cb.set_label("Within-compound $z$-score", fontsize=15)
    cb.ax.tick_params(labelsize=13)
    for ax, L in zip(axes, "AB"):
        ax.text(-0.10, 1.10, L, transform=ax.transAxes, **PANEL)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "figAT3_selectivity.png"))
    plt.close(fig)
    print("figAT3_selectivity.png")


# ----------------------------------------------------------------------
# Figure 4 — concentration-dependent activation (panels A and B only)
# ----------------------------------------------------------------------
def fpl(x, bottom, top, ec50, hill):
    return bottom + (top - bottom) / (1 + (ec50 / x) ** hill)

def fig4_transactivation():
    """Dose-response and 100 uM response, with the two assays drawn apart.

    The PXR and PPAR-alpha assays use different cell backgrounds and reporter
    architectures, so their response magnitudes are not comparable. Earlier
    versions of this figure overlaid them on a shared axis, which invites
    exactly the cross-assay magnitude comparison the manuscript withdraws.
    Each assay is now drawn in its own block with its own axis, and no element
    of the figure places a PXR value next to a PPAR-alpha value.
    """
    d = pd.read_csv(f"{BASE}/data/experimental/transactivation_source_data.csv")
    d["Bio"] = ((d.Replicate - 1) // 3) + 1
    bio = d.groupby(["Compound", "Receptor", "Concentration_uM", "Bio"]).Fold_induction.mean().reset_index()

    COLR = {"PXR": "#1f77b4", "PPAR-a": "#d95f02"}
    LABR = {"PXR": "PXR \u00b7 HepG2 \u00b7 CYP3A4-XREM reporter",
            "PPAR-a": "PPAR-$\\alpha$ \u00b7 COS-1 \u00b7 GAL4-LBD reporter"}

    fig = plt.figure(figsize=(16, 11))
    gs = fig.add_gridspec(3, 4, height_ratios=[1, 1, 1.15],
                          hspace=0.70, wspace=0.34, top=0.93, bottom=0.10,
                          left=0.07, right=0.98)

    # ---- panel A: one ROW per assay, so the two are never on a shared axis ----
    for r, rec in enumerate(["PXR", "PPAR-a"]):
        for c, cmp_ in enumerate(FOUR):
            ax = fig.add_subplot(gs[r, c])
            g = bio[(bio.Compound == cmp_) & (bio.Receptor == rec)]
            m = g.groupby("Concentration_uM").Fold_induction.agg(["mean", "std"])
            ax.errorbar(m.index, m["mean"], yerr=m["std"], fmt="o", ms=7,
                        color=COLR[rec], capsize=4, lw=1.5, zorder=3)
            try:
                p_, _ = curve_fit(fpl, m.index.values.astype(float), m["mean"].values,
                                  p0=[1, m["mean"].max() * 1.2, 20, 1.5], maxfev=200000)
                xs = np.logspace(np.log10(0.3), np.log10(100), 200)
                ax.plot(xs, fpl(xs, *p_), color=COLR[rec], lw=2.1, zorder=2)
            except Exception:
                pass
            ax.set_xscale("log")
            ax.axhline(1, ls=":", color="#999999", lw=1.2)
            # each assay row carries its own scale
            top = max(2.2, m["mean"].max() * 1.35)
            ax.set_ylim(0, top)
            ax.set_title(cmp_, fontsize=16, fontweight="bold", pad=6)
            ax.set_xlabel("Concentration (\u00b5M)", fontsize=14)
            if c == 0:
                ax.set_ylabel("Fold induction", fontsize=15)
            ax.tick_params(labelsize=13)
            ax.spines[["top", "right"]].set_visible(False)
            if c == 0 and r == 0:
                ax.text(-0.30, 1.52, "A", transform=ax.transAxes, **PANEL)

    # row headers naming the assay, one per block
    for y, rec in ((0.965, "PXR"), (0.645, "PPAR-a")):
        fig.text(0.5, y, LABR[rec], ha="center", va="bottom",
                 fontsize=16, fontweight="bold", color=COLR[rec])

    # ---- panel B: two separate blocks, independent axes ----
    hi = bio[bio.Concentration_uM == 100]
    for k, rec in enumerate(["PXR", "PPAR-a"]):
        axB = fig.add_subplot(gs[2, k * 2:k * 2 + 2])
        order = sorted(FOUR, key=lambda c: -hi[(hi.Compound == c) & (hi.Receptor == rec)].Fold_induction.mean())
        mu = [hi[(hi.Compound == c) & (hi.Receptor == rec)].Fold_induction.mean() for c in order]
        sd = [hi[(hi.Compound == c) & (hi.Receptor == rec)].Fold_induction.std(ddof=1) for c in order]
        x = np.arange(len(order))
        axB.bar(x, mu, 0.6, yerr=sd, capsize=5, color=COLR[rec], error_kw=dict(lw=1.6))
        axB.set_xticks(x); axB.set_xticklabels(order, fontsize=15)
        axB.set_ylim(0, max(mu) * 1.30)
        axB.set_ylabel("Fold induction at 100 \u00b5M", fontsize=15)
        axB.set_title(LABR[rec], fontsize=15, fontweight="bold",
                      color=COLR[rec], pad=10)
        axB.tick_params(labelsize=13)
        axB.spines[["top", "right"]].set_visible(False)
        if k == 0:
            axB.text(-0.17, 1.32, "B", transform=axB.transAxes, **PANEL)

    fig.text(0.5, 0.015,
             "The two assays use different cell backgrounds and reporter architectures. "
             "Each is drawn on its own axis; response magnitudes are not comparable between them, "
             "and no comparison between the two is made anywhere in this figure.",
             ha="center", fontsize=14, style="italic", color="#333333")
    fig.subplots_adjust(bottom=0.09)
    fig.savefig(os.path.join(OUT, "figAT4_transactivation.png"))
    plt.close(fig)
    print("figAT4_transactivation.png")


# ----------------------------------------------------------------------
# Supplementary figure — NHANES exposure anchoring
# ----------------------------------------------------------------------
# Replaces the earlier version of this figure, which plotted the withdrawn
# hazard_H literature index on the vertical axis of panel A. hazard_H was not
# auditable (per-component values were never retained) and has been replaced
# throughout by the documentary regulatory priority score, so the figure is
# regenerated here from the deposited data and is reproducible from a clean
# clone like every other figure.
NHANES_LOD = {"PFOS": 99.9, "PFOA": 99.9, "PFHxS": 99.0,
              "PFNA": 95.0, "PFDA": 81.0, "PFUnDA": 59.0}

def figS_nhanes():
    d = pd.read_csv(f"{BASE}/data/docking/analysis_input.csv")
    d["cls"] = d.pfas.map(_klass)
    sub = d[d.pfas.isin(NHANES_LOD)].copy()
    sub["lod"] = sub.pfas.map(NHANES_LOD)
    sub = sub.sort_values("s_bind")
    r6, p6 = stats.spearmanr(sub.s_bind, sub.reg_score)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(14, 6))

    rng = np.random.default_rng(0)
    other = d[~d.pfas.isin(NHANES_LOD)]
    axA.scatter(other.s_bind, other.reg_score + rng.uniform(-.12, .12, len(other)),
                s=70, facecolors="none", edgecolors="#999999", linewidths=1.2,
                label="not an NHANES serum analyte")
    axA.scatter(sub.s_bind, sub.reg_score, s=150, color="#c0392b", zorder=3,
                label="NHANES serum analyte")
    # fixed label offsets: the six analytes are tightly clustered at score 3
    OFF = {"PFHxS": (-10, 14), "PFOA": (-6, -26), "PFOS": (-4, 15),
           "PFDA": (6, -26), "PFUnDA": (10, 14), "PFNA": (10, 8)}
    for _, r in sub.iterrows():
        dx, dy = OFF[r.pfas]
        axA.annotate(r.pfas, (r.s_bind, r.reg_score), textcoords="offset points",
                     xytext=(dx, dy), fontsize=13, color="#c0392b",
                     ha="center", zorder=4)
    axA.set_xlabel("Cross-target binding strength $s_{\\rm bind}$ (kcal mol$^{-1}$)")
    axA.set_ylabel("Regulatory priority score\n(binding instruments naming the compound)")
    axA.set_yticks(range(5))
    axA.legend(frameon=False, fontsize=13, loc="upper left")
    axA.spines[["top", "right"]].set_visible(False)
    axA.text(-0.10, 1.06, "A", transform=axA.transAxes, **PANEL)
    axA.set_ylim(-0.6, 5.0)
    axA.text(0.98, 0.52, f"within the six: $\\rho$ = {r6:+.2f}, p = {p6:.2f}, n = 6",
             transform=axA.transAxes, ha="right", va="top", fontsize=13,
             color="#333333")

    y = np.arange(len(sub))
    axB.barh(y, sub.s_bind, color="#c0392b", alpha=.85, height=.62)
    axB.set_yticks(y); axB.set_yticklabels(sub.pfas)
    axB.set_xlim(0, max(sub.s_bind) * 1.22)
    axB.set_xlabel("Mean |$\\Delta G$| across the seven targets (kcal mol$^{-1}$)")
    for yi, (_, r) in zip(y, sub.iterrows()):
        axB.text(r.s_bind + .12, yi, f"{r.lod:.0f}% > LOD", va="center", fontsize=13)
    axB.spines[["top", "right"]].set_visible(False)
    axB.text(-0.10, 1.06, "B", transform=axB.transAxes, **PANEL)

    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "figATS_nhanes.png"))
    plt.close(fig)
    print(f"figATS_nhanes.png   (n=6 rho={r6:+.3f}, p={p6:.3f})")


if __name__ == "__main__":
    fig1_validation()
    fig2_concordance()
    fig3_selectivity()
    fig4_transactivation()
    figS_nhanes()
    print("\nwritten to", OUT)
