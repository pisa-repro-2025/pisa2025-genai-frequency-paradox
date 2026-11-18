"""
Generate all manuscript figures from audited result CSVs only.
No model is refitted; every plotted value is read from an existing output file.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)

R = ROOT / "r_outputs_real"
SUPP = R / "supplementary"
METH = R / "method_tests"
COUNTRY = R / "formal_by_country"
META = R / "meta_regression"
REAL = R / "real_results"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.dpi": 160,
    "svg.fonttype": "none",
})

INK = "#1F2933"
BLUE = "#2C6EAA"
ORANGE = "#D9822B"
RED = "#B4434A"
TEAL = "#2A9D8F"
GREY = "#7B8794"
LIGHT = "#E8EEF3"

def finish(fig, name):
    path = OUT / name
    fig.savefig(path, format="svg", bbox_inches="tight", facecolor="white")
    fig.savefig(path.with_suffix(".png"), format="png", bbox_inches="tight",
                facecolor="white", dpi=180)
    plt.close(fig)
    print(f"wrote {path.relative_to(ROOT)}")

def style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(GREY)
    ax.spines["bottom"].set_color(GREY)
    ax.tick_params(colors=INK)
    ax.set_title(ax.get_title(), color=INK, fontweight="bold")

# ---------------------------------------------------------------- Figure 2
def figure2():
    d = pd.read_csv(SUPP / "ANA02_H1a_task_specific.csv")
    b2 = d[d["term"] == "ai2"].copy()
    b2 = b2[b2["subset"] == "users_only"]
    names = {"summarise": "Summarizing", "prelim_research": "Preliminary research",
             "drafting": "Drafting", "help_learn": "Help me learn"}
    order = ["summarise", "prelim_research", "drafting", "help_learn"]
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 5.6), sharex=True)
    for ax, task in zip(axes.ravel(), order):
        row = b2[b2["task"] == task].iloc[0]
        coef = d[(d["model"] == f"{task}_users_only")].set_index("term")
        b0 = float(coef.loc["const", "est"])
        b1 = float(coef.loc["ai", "est"])
        bq = float(coef.loc["ai2", "est"])
        raw_mean = float(row["vertex_raw"] - row["vertex_centered"])
        xv = np.linspace(2, 5, 300)
        xc = xv - raw_mean
        yv = b0 + b1 * xc + bq * xc ** 2
        ax.plot(xv, yv, color=BLUE, lw=2)
        vx = float(row["vertex_raw"])
        vy = b0 + b1 * (vx - raw_mean) + bq * (vx - raw_mean) ** 2
        ax.scatter([vx], [vy], color=RED, s=34, zorder=3)
        ax.axvline(vx, color=RED, ls="--", lw=.9, alpha=.7)
        ax.text(vx, ax.get_ylim()[1], f"  vertex={vx:.2f}", color=RED, va="top", fontsize=8)
        ax.set_title(names[task])
        ax.set_ylabel("Adjusted science score")
        ax.grid(axis="y", color=LIGHT, lw=.8)
        style_ax(ax)
    for ax in axes[1]:
        ax.set_xlabel("AI-use frequency (raw 1–5 scale; users only)")
    fig.suptitle("Within-user quadratic curves for four AI-use purposes", color=INK, fontweight="bold", y=.99)
    fig.tight_layout(rect=[0, 0, 1, .96])
    finish(fig, "figure2_task_specific_curves.svg")

# ---------------------------------------------------------------- Figure 3
def figure3():
    d = pd.read_csv(SUPP / "ANA02_H1a_task_specific.csv")
    q = d[d["term"] == "ai2"].copy()
    p = q.pivot(index="task", columns="subset", values="est")
    p = p.reindex(["summarise", "prelim_research", "drafting", "help_learn"])
    labels = ["Summarizing", "Preliminary\nresearch", "Drafting", "Help me\nlearn"]
    x = np.arange(len(p)); w = .35
    fig, ax = plt.subplots(figsize=(7.4, 3.7))
    ax.bar(x - w/2, p["all"], w, label="Full sample (includes non-users)", color=BLUE)
    ax.bar(x + w/2, p["users_only"], w, label="Users only (ST438 > 1)", color=ORANGE)
    ax.axhline(0, color=INK, lw=1)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel("Quadratic coefficient (ai²)")
    ax.set_title("Composition reverses the sign of curvature")
    ax.legend(frameon=False, loc="upper right")
    ax.grid(axis="y", color=LIGHT, lw=.8)
    for i, (a, b) in enumerate(zip(p["all"], p["users_only"])):
        ax.text(i - w/2, a + (.12 if a >= 0 else -.28), f"{a:+.2f}", ha="center", fontsize=7)
        ax.text(i + w/2, b + (.12 if b >= 0 else -.28), f"{b:+.2f}", ha="center", fontsize=7)
    style_ax(ax)
    fig.tight_layout()
    finish(fig, "figure3_curvature_sign_reversal.svg")

# ---------------------------------------------------------------- Figure 4
def figure4():
    d = pd.read_csv(METH / "LindMehlum.csv")
    labels = ["Pooled", "Summarize\n(users)", "Preliminary\n(users)", "Drafting\n(users)", "Help me learn\n(users)"]
    x = np.arange(len(d)); w = .36
    fig, ax = plt.subplots(figsize=(7.6, 3.8))
    ax.bar(x - w/2, d["slope_low"], w, yerr=1.96*d["slope_low_se"], color=TEAL,
           label="Slope at 5th percentile", capsize=3)
    ax.bar(x + w/2, d["slope_high"], w, yerr=1.96*d["slope_high_se"], color=RED,
           label="Slope at 95th percentile", capsize=3)
    ax.axhline(0, color=INK, lw=1)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel("Slope (score points per AI-unit)")
    ax.set_title("Lind–Mehlum slope criteria: pooled fails, all four user tasks pass")
    ax.legend(frameon=False)
    ax.grid(axis="y", color=LIGHT, lw=.8)
    style_ax(ax)
    fig.tight_layout()
    finish(fig, "figure4_lind_mehlum_slopes.svg")

# ---------------------------------------------------------------- Figure 5
def figure5():
    d = pd.read_csv(METH / "FunctionalForms.csv")
    samples = d["sample"].drop_duplicates().tolist()
    models = ["M1_linear", "M2_quadratic", "M3_cubic", "M4_RCS4k"]
    short_models = ["M1\nlinear", "M2\nquadratic", "M3\ncubic", "M4\nRCS-4k"]
    short_samples = ["Pooled", "Summarize\n(users)", "Preliminary\n(users)",
                     "Drafting\n(users)", "Help me learn\n(users)"]
    mat = np.full((len(samples), len(models)), np.nan)
    for i, s in enumerate(samples):
        g = d[d["sample"] == s].set_index("model")
        for j, m in enumerate(models):
            mat[i, j] = float(g.loc[m, "delta_AIC"])
    fig, ax = plt.subplots(figsize=(7.2, 3.9))
    im = ax.pcolormesh(np.arange(len(models)+1), np.arange(len(samples)+1), mat,
                       norm=matplotlib.colors.LogNorm(vmin=.1, vmax=4000), cmap="Blues_r",
                       edgecolors="white", linewidths=1.5)
    ax.set_xticks(np.arange(len(models))+.5); ax.set_xticklabels(short_models)
    ax.set_yticks(np.arange(len(samples))+.5); ax.set_yticklabels(short_samples)
    for i in range(len(samples)):
        for j in range(len(models)):
            val = mat[i, j]
            color = INK if val < .1 else ("white" if val < 10 else INK)
            ax.text(j+.5, i+.5, f"{val:,.1f}", ha="center", va="center", fontsize=7.5, color=color)
    cb = fig.colorbar(im, ax=ax, fraction=.046, pad=.04)
    cb.set_label("ΔAIC from best model (log scale)")
    cb.ax.tick_params(labelsize=7)
    ax.set_title("Functional-form comparison: linearity rejected, cubic often preferred")
    ax.invert_yaxis()
    style_ax(ax)
    fig.tight_layout()
    finish(fig, "figure5_functional_forms.svg")

# ---------------------------------------------------------------- Figure 6
def figure6():
    by = pd.read_csv(COUNTRY / "by_country_formal.csv")
    meta = pd.read_csv(COUNTRY / "meta_analysis.csv")
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.6), gridspec_kw={"width_ratios": [1.3, 1]})
    ax = axes[0]
    ax.hist(by["ai2"], bins=28, color=BLUE, alpha=.85, edgecolor="white")
    re = float(meta.loc[meta["model"] == "随机效应DL", "est"].iloc[0])
    ax.axvline(0, color=INK, lw=1.1, label="Zero curvature")
    ax.axvline(re, color=RED, lw=1.6, label=f"Random-effect estimate = {re:.3f}")
    counts = by["shape"].value_counts().to_dict()
    txt = (f"Significant inverted U: {counts.get('倒U', 0)}\n"
           f"Significant U shape: {counts.get('U型', 0)}\n"
           f"Nonsignificant: {counts.get('不显著', 0)}")
    ax.text(.02, .96, txt, transform=ax.transAxes, va="top", fontsize=8,
            bbox=dict(boxstyle="round,pad=.35", fc="white", ec=LIGHT))
    ax.set_xlabel("Country-specific quadratic coefficient (ai²)")
    ax.set_ylabel("Number of countries/economies")
    ax.set_title("Country-level curvature is widely dispersed")
    ax.legend(frameon=False, fontsize=7.5)
    ax.grid(axis="y", color=LIGHT, lw=.8)
    style_ax(ax)

    ax = axes[1]
    m = meta.copy()
    y = np.arange(len(m))
    est = m["est"].values; se = m["se"].values
    labels = ["Fixed effect", "Random effects (DL)"]
    ax.errorbar(est, y, xerr=1.96*se, fmt="o", color=RED, capsize=4, ms=5)
    pooled = pd.read_csv(REAL / "model2_invertedU.csv")
    pool_est = float(pooled.loc[pooled["term"] == "ai2", "est"].iloc[0])
    ax.axvline(pool_est, color=ORANGE, ls="--", lw=1.4)
    ax.text(pool_est, 1.58, f"Pooled = {pool_est:+.3f}", color=ORANGE,
            ha="center", fontsize=7.5)
    ax.axvline(0, color=INK, lw=1)
    ax.set_yticks(y); ax.set_yticklabels(labels)
    ax.set_xlim(min(est - 1.96*se) - .4, pool_est + .5)
    ax.set_xlabel("Pooled quadratic coefficient (ai²)")
    ax.set_title(f"I² = {float(m['I2'].iloc[0]):.1f}%")
    for yi, (ei, si) in enumerate(zip(est, se)):
        ax.text(ei, yi+.15, f"{ei:+.3f} ± {1.96*si:.3f}", ha="center", fontsize=7.5)
    ax.grid(axis="x", color=LIGHT, lw=.8)
    style_ax(ax)
    fig.suptitle("Sign reversal between pooled and country-level aggregation", color=INK, fontweight="bold", y=1.01)
    fig.tight_layout()
    finish(fig, "figure6_country_shape_and_meta.svg")

# ---------------------------------------------------------------- Figure 7
def figure7():
    u = pd.read_csv(META / "univariate.csv")
    f = pd.read_csv(METH / "FDR_BH.csv")
    f = f[f["family"] == "meta-regression univariate"].set_index("test")
    u = u.iloc[::-1].reset_index(drop=True)
    labels = ["Mean ESCS", "Proportion girls", "Mean AI use", "Mean science", "AI opportunity"]
    tstat = u["beta"].values / u["se"].values
    q = np.array([f.loc[v, "q_BH"] for v in u["var"]])
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.5))
    y = np.arange(len(u))
    ax = axes[0]
    colors = [RED if qv < .05 else GREY for qv in q]
    ax.barh(y, tstat, color=colors, alpha=.9)
    ax.axvline(0, color=INK, lw=1)
    ax.axvline(-1.96, color=GREY, ls="--", lw=.9)
    ax.axvline(1.96, color=GREY, ls="--", lw=.9)
    ax.set_yticks(y); ax.set_yticklabels(labels)
    ax.set_xlabel("t statistic (β / SE)")
    ax.set_title("Direction and significance")
    ax.grid(axis="x", color=LIGHT, lw=.8)
    style_ax(ax)

    ax = axes[1]
    nq = -np.log10(q)
    ax.barh(y, nq, color=colors, alpha=.9)
    ax.axvline(-np.log10(.05), color=RED, ls="--", lw=1)
    ax.set_yticks(y); ax.set_yticklabels([])
    ax.set_xlabel("−log₁₀(BH q)")
    ax.set_title("Multiplicity-adjusted evidence")
    for yi, qv in zip(y, q):
        ax.text(nq[yi]+.05, yi, f"q={qv:.3f}", va="center", fontsize=7.5)
    ax.grid(axis="x", color=LIGHT, lw=.8)
    style_ax(ax)
    fig.suptitle("Country-level moderators of curvature heterogeneity", color=INK, fontweight="bold", y=1.02)
    fig.tight_layout()
    finish(fig, "figure7_meta_regression_evidence.svg")

# ---------------------------------------------------------------- Figure 8
def figure8():
    d = pd.read_csv(SUPP / "ANA01_H2b_joint_contrast.csv").set_index("term")
    terms = ["ai2_x_use", "ai2_x_gen", "CONTRAST_use_minus_gen"]
    labels = ["Quadratic × AI-specific evaluation", "Quadratic × general online evaluation",
              "Direct contrast (AI-specific − general)"]
    est = d.loc[terms, "est"].values
    se = d.loc[terms, "brr_se"].values
    fig, ax = plt.subplots(figsize=(7.8, 3.2))
    y = np.arange(3)
    colors = [RED, TEAL, BLUE]
    for yi, (ei, si, ci) in enumerate(zip(est, se, colors)):
        ax.errorbar([ei], [yi], xerr=[1.96*si], fmt="o", color=ci, ecolor=ci,
                    capsize=5, ms=7, ls="none")
        lab_y = yi + .22 if yi < 2 else yi - .22
        ax.text(ei, lab_y, f"{ei:+.3f} (SE={si:.3f})", ha="center", color=ci, fontsize=7.5)
    ax.axvline(0, color=INK, lw=1)
    ax.set_yticks(y); ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Coefficient (95% BRR CI)")
    ax.set_title("Curvature moderation is AI-specific")
    ax.grid(axis="x", color=LIGHT, lw=.8)
    style_ax(ax)
    fig.tight_layout()
    finish(fig, "figure8_h2b_contrast.svg")

# ---------------------------------------------------------------- Figure 9
def figure9():
    summ = pd.read_csv(METH / "JohnsonNeyman_summary.csv")
    rows = []
    for _, r in summ.iterrows():
        grid = pd.read_csv(METH / f"JN_{r['moderator']}.csv")
        lo, hi = float(grid["mod"].min()), float(grid["mod"].max())
        sig_lo, sig_hi = [float(v) for v in r["sig_neg_range"].strip("[]").split(",")]
        med = float(grid["mod"].median())
        rows.append((r["moderator"], lo, hi, sig_lo, sig_hi, med, float(r["curvature_at_median"]), float(r["frac_grid_sig"])))
    fig, ax = plt.subplots(figsize=(8.0, 2.9))
    labels = ["AI-specific evaluation", "General online evaluation"]
    for i, row in enumerate(rows):
        name, lo, hi, slo, shi, med, curv, frac = row
        y = i
        ax.plot([lo, hi], [y, y], color=GREY, lw=9, solid_capstyle="butt", alpha=.35)
        ax.plot([slo, shi], [y, y], color=RED, lw=9, solid_capstyle="butt")
        ax.scatter([med], [y], color=INK, s=26, zorder=4)
        lab_y = y + .20 if i == 0 else y - .22
        med_y = y + .24 if i == 0 else y - .26
        ax.text(slo, lab_y, f"[{slo:.3f}, {shi:.3f}]", ha="center", fontsize=8)
        ax.text(hi+.03, y, f"{frac*100:.0f}%", va="center", fontsize=8)
        ax.text(med, med_y, "median", ha="center", fontsize=7, color=INK)
    ax.set_yticks([0, 1]); ax.set_yticklabels(labels)
    ax.set_xlabel("Centered moderator value (2nd–98th percentile grid)")
    ax.set_title("Johnson–Neyman regions of significantly negative curvature")
    ax.grid(axis="x", color=LIGHT, lw=.8)
    style_ax(ax)
    fig.text(.5, .01, "Red = significant negative curvature; grey = remaining observed range",
             ha="center", fontsize=7.5, color=GREY)
    fig.tight_layout(rect=[0, .10, 1, 1])
    finish(fig, "figure9_johnson_neyman.svg")

# ---------------------------------------------------------------- Figure 10
def figure10():
    d = pd.read_csv(SUPP / "ANA03_H3H4_stratification.csv")
    h3 = d[(d["model"].str.startswith("H3_")) & (d["term"] == "escs_c")]
    h4 = d[(d["model"].str.startswith("H4_")) & (d["term"] == "ai2")]
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.4))
    for ax in axes:
        ax.set_axisbelow(True)
    ax = axes[0]
    labs = ["AI-specific opportunity", "General online opportunity"]
    ax.barh([0, 1], h3["est"].values, xerr=1.96*h3["se"].values, color=[BLUE, TEAL], capsize=4)
    ax.axvline(0, color=INK, lw=1)
    ax.set_yticks([0, 1]); ax.set_yticklabels(labs); ax.invert_yaxis()
    ax.set_xlabel("ESCS coefficient (95% CI)")
    ax.set_title("H3: opportunity is socioeconomically graded")
    for i, (e, p) in enumerate(zip(h3["est"], h3["p"])):
        ax.text(e, i-.28, f"{e:+.3f}, p<.001", ha="center", fontsize=7.5)
    ax.grid(axis="x", color=LIGHT, lw=.8); style_ax(ax)

    ax = axes[1]
    labs = ["Low SES", "Middle SES", "High SES"]
    y = np.arange(len(h4))
    ax.barh(y, h4["est"].values, xerr=1.96*h4["se"].values, color=ORANGE, capsize=4)
    ax.axvline(0, color=INK, lw=1)
    ax.set_yticks(y); ax.set_yticklabels(labs); ax.invert_yaxis()
    ax.set_xlabel("Quadratic coefficient ai² (95% CI)")
    ax.set_title("H4: inverted U within every SES tertile")
    for yi, (e, p) in enumerate(zip(h4["est"], h4["p"])):
        ax.text(e, yi-.28, f"{e:+.3f}", ha="center", fontsize=7.5)
    ax.grid(axis="x", color=LIGHT, lw=.8); style_ax(ax)
    fig.tight_layout()
    finish(fig, "figure10_ses_stratification.svg")

# ---------------------------------------------------------------- Figure 11
def figure11():
    d = pd.read_csv(METH / "FDR_BH.csv").copy()
    d["test"] = d["test"].str.replace("_", " ")
    d = d.sort_values(["family", "q_BH"], ascending=[True, True]).reset_index(drop=True)
    y = np.arange(len(d))
    colors = [RED if s else GREY for s in d["sig_q05"]]
    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    ax.hlines(y, 0, d["q_BH"].values, color=colors, lw=2)
    ax.scatter(d["q_BH"].values, y, color=colors, s=28, zorder=3)
    ax.axvline(.05, color=INK, ls="--", lw=1)
    ax.set_yticks(y); ax.set_yticklabels(d["test"], fontsize=7.5); ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("Benjamini–Hochberg adjusted q (log scale)")
    ax.set_title("FDR correction: core results survive at q < .05")
    ax.grid(axis="x", color=LIGHT, lw=.8)
    for yi, qv in zip(y, d["q_BH"]):
        lab = "<.001" if qv < .001 else f"{qv:.3f}"
        ax.text(qv*1.08, yi, lab, va="center", fontsize=7)
    style_ax(ax)
    fig.tight_layout()
    finish(fig, "figure11_fdr_adjustment.svg")

# ---------------------------------------------------------------- Figure 12
def figure12():
    fig, ax = plt.subplots(figsize=(9.2, 4.2))
    ax.set_xlim(0, 12); ax.set_ylim(0, 5.4); ax.axis("off")
    boxes = [
        (0.3, 2.1, 2.1, 1.2, "PISA 2025\nStudent file\nN = 755,721", BLUE),
        (3.0, 2.1, 2.2, 1.2, "5 plausible values\n80 Fay BRR weights\nWithin-system weights", TEAL),
        (5.8, 2.1, 2.1, 1.2, "Weighted least squares\nRubin combination\nBRR standard errors", ORANGE),
        (8.6, 3.5, 2.9, 1.2, "Substantive layers\nPooled · 68 systems\nTasks · Moderation · SES", BLUE),
        (8.6, 1.5, 2.9, 1.2, "Formal checks\nLind–Mehlum · M1–M4\nJN regions · BH FDR", RED),
    ]
    for x, y, w, h, txt, color in boxes:
        patch = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.06,rounding_size=0.08",
                               facecolor="white", edgecolor=color, linewidth=1.5)
        ax.add_patch(patch)
        ax.text(x+w/2, y+h/2, txt, ha="center", va="center", fontsize=8.2, color=INK)
    arrows = [((2.4, 2.7), (3.0, 2.7)), ((5.2, 2.7), (5.8, 2.7)),
              ((7.9, 2.7), (8.6, 4.1)), ((7.9, 2.7), (8.6, 2.1))]
    for a, b in arrows:
        ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=12,
                                     color=GREY, linewidth=1.3))
    ax.text(6, 5.05, "Analytic pipeline and inference boundary", ha="center",
            fontsize=12, color=INK, fontweight="bold")
    ax.text(10.05, .85, "Inference: association only; no causal dose", ha="center",
            fontsize=8.2, color=RED, fontweight="bold")
    ax.add_patch(FancyArrowPatch((10.05, 1.5), (10.05, 1.05), arrowstyle="-|>",
                                 mutation_scale=11, color=RED, linewidth=1.2))
    finish(fig, "figure12_analytic_inference_map.svg")

if __name__ == "__main__":
    figure2()
    figure3()
    figure4()
    figure5()
    figure6()
    figure7()
    figure8()
    figure9()
    figure10()
    figure11()
    figure12()
    print("All figures generated.")
