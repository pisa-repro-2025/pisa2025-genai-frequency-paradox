# -*- coding: utf-8 -*-
"""
方法强化批次 1：四项正式检验，全部沿用主分析同一套 WLS + Rubin + Fay BRR 管道
  (A) Lind–Mehlum (2010) 倒U正式检验：低/高端斜率 + 顶点位置 + BRR SE
  (B) 功能形式比较 M1–M4：线性 / 二次 / 三次 / 限制立方样条(RCS)，加权 AIC/BIC
  (C) Johnson–Neyman 显著性区域：H2b 联合模型中调节变量在何处使曲率为负/显著
  (D) FDR (Benjamini–Hochberg) 多重比较校正：四任务×两样本 + 元回归单变量

立场：如实报告，不预设、不迎合论文叙事。scipy 的 optimize/stats 可加载（statsmodels 被拦），
      故检验用纯 numpy + scipy.stats 实现，数学等价。
"""
import pyreadstat, numpy as np, pandas as pd
from scipy import stats
import os, time, json
t0 = time.time()
OUT = "r_outputs_real/method_tests"
os.makedirs(OUT, exist_ok=True)

pv = [f"PV{i}SCIE" for i in range(1, 6)]
rw = [f"W_FSTURWT{i}" for i in range(1, 81)]
rho = 0.5
R = len(rw)

# ---------------------------------------------------------------- 读数据
cols = ["CNT", "AIUSESCH", "ST438Q01DA", "ST438Q02DA", "ST438Q03DA", "ST438Q04DA",
        "ST436Q16DA", "ST436Q02DA", "ESCS", "ST004D01T"] + pv + ["W_FSTUWT"] + rw
print("reading sav...", flush=True)
df, _ = pyreadstat.read_sav("data/CY09_MS_STU_PUF.sav", usecols=cols, apply_value_formats=False)
print(f"N={len(df):,} read_sec={time.time()-t0:.1f}", flush=True)

for c in ["ST436Q16DA", "ST436Q02DA"]:
    df.loc[df[c].isin([5, 6, 7, 8, 9]), c] = np.nan
for c in ["ST438Q01DA", "ST438Q02DA", "ST438Q03DA", "ST438Q04DA"]:
    df.loc[df[c].isin([95, 96, 97, 98, 99]), c] = np.nan
df.loc[df["ESCS"].isin([95, 96, 97, 98, 99]), "ESCS"] = np.nan
df.loc[~df["ST004D01T"].isin([1, 2]), "ST004D01T"] = np.nan
df.loc[(~np.isfinite(df["AIUSESCH"])) | (df["AIUSESCH"].abs() > 10), "AIUSESCH"] = np.nan
for c in pv:
    df.loc[df[c] >= 9997, c] = np.nan
df["gender"] = df["ST004D01T"] - 1

# ---------------------------------------------------------------- 工具函数
def norm_weights(w, countries):
    s = pd.Series(w); g = pd.Series(countries)
    return w * s.groupby(g).transform(lambda x: len(x)/x.sum()).values

def wls(X, y, w):
    sw = np.sqrt(w); Xw = X*sw[:, None]; yw = y*sw
    beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    return beta

def wls_cov(X, y, w):
    beta = wls(X, y, w)
    resid = y - X@beta
    Z = X*w[:, None]*resid[:, None]
    meat = Z.T@Z
    inv = np.linalg.pinv(X.T@(X*w[:, None]))
    return beta, inv@meat@inv

def brr_coefs(X, y, w_full, rwm):
    """返回 (R, k) 的 BRR 重复系数矩阵"""
    bs = np.empty((R, X.shape[1]))
    for r in range(R):
        wr = (1-rho)*w_full + rho*rwm[:, r]
        bs[r] = wls(X, y, wr)
    return bs

def brr_var_from_draws(draws):
    """draws: (R,) 或 (R,k) -> BRR 方差（标量或向量）"""
    draws = np.asarray(draws)
    return ((draws - draws.mean(0))**2).sum(0)/(R*(1-rho)**2)

def combine_across_pv(point, var_within):
    """Rubin 合并：point (M,) 每个PV的点估计；var_within (M,) 每个PV的BRR方差
       返回 qbar, se, t, p, df"""
    point = np.asarray(point, float); var_within = np.asarray(var_within, float)
    M = len(point); qbar = point.mean(); ubar = var_within.mean()
    b = ((point-qbar)**2).sum()/(M-1) if M > 1 else 0.0
    T = ubar + (1+1/M)*b
    if b > 0 and np.isfinite(b):
        dfree = (M-1)*(1+ubar/((1+1/M)*b))**2
    else:
        dfree = 1e6
    se = np.sqrt(max(T, 0))
    t = qbar/se if se > 0 else np.nan
    p = 2*stats.t.sf(abs(t), dfree) if np.isfinite(t) else np.nan
    return qbar, se, t, p, dfree

def prep(dsub):
    """返回 w, rwm（国家内归一化权重）"""
    w = norm_weights(dsub["W_FSTUWT"].values, dsub["CNT"].values)
    rwm = np.column_stack([norm_weights(dsub[r].values, dsub["CNT"].values) for r in rw])
    return w, rwm

def make_sub(need):
    return df.dropna(subset=need+pv+["W_FSTUWT"]+rw).copy()

# ================================================================ (A) Lind–Mehlum
def lind_mehlum(dsub, x, label, raw_mean=None, lo_q=0.05, hi_q=0.95):
    """x: 已中心化的 AI 变量列名（设计中的 ai）。检验低/高端斜率与顶点。
       顶点与斜率均用 BRR 重复系数 delta 法（直接在每个重复中计算），Rubin 跨 PV 合并。
       raw_mean: 该子样本中原始题(未中心化)的均值，用于把顶点还原到 1–5 原始量表。"""
    d = dsub.copy()
    d["x"] = d[x].values; d["x2"] = d["x"]**2
    d["escs_c"] = d["ESCS"] - d["ESCS"].mean()
    X = pd.DataFrame({"const": 1.0, "x": d["x"].values, "x2": d["x2"].values,
                      "escs_c": d["escs_c"].values, "gender": d["gender"].values})
    Xm = X.values; names = list(X.columns)
    ix1, ix2 = names.index("x"), names.index("x2")
    w, rwm = prep(d)
    xlo = np.quantile(d["x"].values, lo_q); xhi = np.quantile(d["x"].values, hi_q)
    xmin, xmax = d["x"].min(), d["x"].max()
    # 每个 PV：点估计 + BRR 重复系数
    draws_all = []; points_all = []
    for pvc in pv:
        y = d[pvc].values.astype(float)
        b = wls(Xm, y, w)
        bs = brr_coefs(Xm, y, w, rwm)
        draws_all.append(bs); points_all.append(b)
    draws_all = np.array(draws_all)   # (M, R, k)
    points_all = np.array(points_all) # (M, k)
    # 斜率函数 s(t) = b1 + 2*b2*t；顶点 = -b1/(2*b2)
    def slope_stats(t):
        pts = points_all[:, ix1] + 2*points_all[:, ix2]*t
        dr = draws_all[:, :, ix1] + 2*draws_all[:, :, ix2]*t   # (M,R)
        varw = np.array([brr_var_from_draws(dr[m]) for m in range(len(pv))])
        return combine_across_pv(pts, varw)
    s_lo = slope_stats(xlo); s_hi = slope_stats(xhi)
    vert_pts = -points_all[:, ix1]/(2*points_all[:, ix2])
    vert_draw = -draws_all[:, :, ix1]/(2*draws_all[:, :, ix2])
    varw_v = np.array([brr_var_from_draws(vert_draw[m]) for m in range(len(pv))])
    v = combine_across_pv(vert_pts, varw_v)
    # 顶点还原到 raw 尺度：用原始题(未中心化)均值，而非中心化 ai 的均值(≈0)
    x_mean = raw_mean if raw_mean is not None else d[x].mean()
    crit_low = (s_lo[0] > 0) and (s_lo[3] < 0.05)   # 低端斜率显著正
    crit_high = (s_hi[0] < 0) and (s_hi[3] < 0.05)  # 高端斜率显著负
    crit_vertex = (v[0] >= xlo) and (v[0] <= xhi)
    verdict = "inverted-U supported" if (crit_low and crit_high and crit_vertex) else "not fully supported"
    return dict(label=label, N=len(d), x_mean=x_mean,
                x_lo_centered=xlo, x_hi_centered=xhi,
                x_lo_raw=xlo+x_mean, x_hi_raw=xhi+x_mean,
                slope_low=s_lo[0], slope_low_se=s_lo[1], slope_low_p=s_lo[3],
                slope_high=s_hi[0], slope_high_se=s_hi[1], slope_high_p=s_hi[3],
                vertex_centered=v[0], vertex_se=v[1], vertex_p=v[3],
                vertex_raw=v[0]+x_mean,
                vertex_in_range=bool(crit_vertex),
                crit_low=bool(crit_low), crit_high=bool(crit_high),
                verdict=verdict)

print("\n[A] Lind–Mehlum 正式倒U检验", flush=True)
lm_rows = []
# A1: pooled 主分析（AIUSESCH）
d_main = make_sub(["AIUSESCH", "ESCS", "ST004D01T"])
d_main["ai"] = d_main["AIUSESCH"] - d_main["AIUSESCH"].mean()
lm_rows.append(lind_mehlum(d_main, "ai", "Pooled (AIUSESCH, full sample)",
                           raw_mean=d_main["AIUSESCH"].mean()))
print(f"  Pooled: {lm_rows[-1]['verdict']}  slope_lo={lm_rows[-1]['slope_low']:+.3f}(p={lm_rows[-1]['slope_low_p']:.2g}) "
      f"slope_hi={lm_rows[-1]['slope_high']:+.3f}(p={lm_rows[-1]['slope_high_p']:.2g}) "
      f"vertex_raw={lm_rows[-1]['vertex_raw']:.2f}", flush=True)
# A2: 四任务用户内（v2 核心 H1a）
task_cols = ["ST438Q01DA", "ST438Q02DA", "ST438Q03DA", "ST438Q04DA"]
task_names = ["summarise", "preliminary_research", "drafting", "help_me_learn"]
for tc, tn in zip(task_cols, task_names):
    d = make_sub([tc, "ESCS", "ST004D01T"])
    d = d[d[tc] > 1].copy()          # 仅用户（与 ANA-02 完全一致）
    raw_mean = d[tc].mean()          # 原始 1–5 题均值，用于还原顶点
    d["ai"] = d[tc] - raw_mean
    r = lind_mehlum(d, "ai", f"Users-only: {tn}", raw_mean=raw_mean)
    lm_rows.append(r)
    print(f"  {tn:20s}: {r['verdict']:20s} slope_lo={r['slope_low']:+.3f}(p={r['slope_low_p']:.1g}) "
          f"slope_hi={r['slope_high']:+.3f}(p={r['slope_high_p']:.1g}) vertex_raw={r['vertex_raw']:.2f}", flush=True)
LM = pd.DataFrame(lm_rows)
LM.to_csv(f"{OUT}/LindMehlum.csv", index=False, encoding="utf-8")

# ================================================================ (B) 功能形式比较
def rcs_basis(x, knots):
    """限制立方样条基（Harrell）：返回 (n, m-1) = [x, nu_1..nu_{m-2}]"""
    k = np.sort(np.asarray(knots, float)); m = len(k)
    t = (x - k[0])/(k[-1]-k[0]); tk = (k-k[0])/(k[-1]-k[0])
    xb = x.copy()
    cols = [xb]
    for j in range(1, m-1):
        nu = (np.clip(t-tk[j], 0, None)**3
              - np.clip(t-tk[m-2], 0, None)**3*(tk[m-1]-tk[j])/(tk[m-1]-tk[m-2])
              + np.clip(t-tk[m-1], 0, None)**3*(tk[m-2]-tk[j])/(tk[m-1]-tk[m-2]))
        cols.append(nu)
    return np.column_stack(cols)

def fit_forms(dsub, x, label):
    d = dsub.copy(); xx = d[x].values
    d["escs_c"] = d["ESCS"] - d["ESCS"].mean()
    base_ctrl = [d["escs_c"].values, d["gender"].values]
    knots = np.quantile(xx, [0.05, 0.35, 0.65, 0.95])
    forms = {
        "M1_linear":   np.column_stack([np.ones(len(d)), xx] + base_ctrl),
        "M2_quadratic":np.column_stack([np.ones(len(d)), xx, xx**2] + base_ctrl),
        "M3_cubic":    np.column_stack([np.ones(len(d)), xx, xx**2, xx**3] + base_ctrl),
        "M4_RCS4k":    np.column_stack([np.ones(len(d)), rcs_basis(xx, knots)] + base_ctrl),
    }
    w, _ = prep(d)
    rows = []
    for nm, X in forms.items():
        # 逐 PV 拟合，取平均加权 RSS 与参数个数
        rss = []
        for pvc in pv:
            y = d[pvc].values.astype(float)
            b = wls(X, y, w); e = y - X@b
            rss.append((w*e**2).sum())
        rssbar = np.mean(rss); n = len(d); k = X.shape[1]
        ll = -0.5*(n*np.log(2*np.pi*rssbar/n) + n - np.log(w).sum())
        aic = -2*ll + 2*k; bic = -2*ll + k*np.log(n)
        rows.append(dict(sample=label, model=nm, k=k, wRSS=rssbar,
                         AIC=aic, BIC=bic, delta_AIC=aic-min(r['AIC'] for r in rows) if rows else 0.0))
    dfm = pd.DataFrame(rows)
    dfm["delta_AIC"] = dfm["AIC"] - dfm["AIC"].min()
    dfm["delta_BIC"] = dfm["BIC"] - dfm["BIC"].min()
    dfm["best_by_AIC"] = dfm["AIC"] == dfm["AIC"].min()
    return dfm

print("\n[B] 功能形式比较 M1–M4（加权 AIC/BIC）", flush=True)
form_rows = []
form_rows.append(fit_forms(d_main, "ai", "Pooled (AIUSESCH)"))
for tc, tn in zip(task_cols, task_names):
    d = make_sub([tc, "ESCS", "ST004D01T"]); d = d[d[tc] > 1].copy()
    d["ai"] = d[tc] - d[tc].mean()
    form_rows.append(fit_forms(d, "ai", f"Users-only: {tn}"))
FORMS = pd.concat(form_rows, ignore_index=True)
FORMS.to_csv(f"{OUT}/FunctionalForms.csv", index=False, encoding="utf-8")
for s, g in FORMS.groupby("sample", sort=False):
    best = g.loc[g.AIC.idxmin()]
    print(f"  {s:26s} best={best['model']:12s} ΔAIC(quad-linear)={g.loc[g.model=='M1_linear','delta_AIC'].iloc[0]:.1f}", flush=True)

# ================================================================ (C) Johnson–Neyman
print("\n[C] Johnson–Neyman 显著性区域（H2b 联合模型）", flush=True)
need_c = ["AIUSESCH", "ST436Q16DA", "ST436Q02DA", "ESCS", "ST004D01T"]
dc = make_sub(need_c)
dc["ai"] = dc["AIUSESCH"] - dc["AIUSESCH"].mean(); dc["ai2"] = dc["ai"]**2
dc["ai_use"] = dc["ST436Q16DA"] - dc["ST436Q16DA"].mean()
dc["ai_gen"] = dc["ST436Q02DA"] - dc["ST436Q02DA"].mean()
dc["escs_c"] = dc["ESCS"] - dc["ESCS"].mean()
X = pd.DataFrame({"const": 1., "ai": dc["ai"].values, "ai2": dc["ai2"].values,
    "ai_use": dc["ai_use"].values, "ai_gen": dc["ai_gen"].values,
    "ai_x_use": dc["ai"].values*dc["ai_use"].values,
    "ai2_x_use": dc["ai2"].values*dc["ai_use"].values,
    "ai_x_gen": dc["ai"].values*dc["ai_gen"].values,
    "ai2_x_gen": dc["ai2"].values*dc["ai_gen"].values,
    "escs_c": dc["escs_c"].values, "gender": dc["gender"].values})
names = list(X.columns); Xm = X.values
i_base = names.index("ai2"); i_use = names.index("ai2_x_use"); i_gen = names.index("ai2_x_gen")
w, rwm = prep(dc)
draws = []; pts = []
for pvc in pv:
    y = dc[pvc].values.astype(float)
    pts.append(wls(Xm, y, w)); draws.append(brr_coefs(Xm, y, w, rwm))
pts = np.array(pts); draws = np.array(draws)

def jn_other(mod_col_idx, other_name):
    """曲率 = b_ai2 + b_interaction * mod；扫描 mod 范围，找曲率显著<0 的区域"""
    lo, hi = np.quantile(dc[other_name].values, [0.02, 0.98])
    grid = np.linspace(lo, hi, 121)
    recs = []
    for g in grid:
        p_pts = pts[:, i_base] + pts[:, mod_col_idx]*g
        p_dr = draws[:, :, i_base] + draws[:, :, mod_col_idx]*g
        varw = np.array([brr_var_from_draws(p_dr[m]) for m in range(len(pv))])
        q, se, t, p, df = combine_across_pv(p_pts, varw)
        recs.append(dict(mod=g, curvature=q, se=se, p=p, sig_neg=bool((q < 0) and (p < 0.05))))
    J = pd.DataFrame(recs)
    J.to_csv(f"{OUT}/JN_{other_name}.csv", index=False, encoding="utf-8")
    sig = J[J.sig_neg]
    if len(sig):
        return dict(moderator=other_name, sig_neg_range=f"[{sig['mod'].min():.3f}, {sig['mod'].max():.3f}]",
                    frac_grid_sig=len(sig)/len(J), curvature_at_median=float(J.loc[(J['mod']-J['mod'].median()).abs().idxmin(), 'curvature']))
    return dict(moderator=other_name, sig_neg_range="none in observed 2–98% range",
                frac_grid_sig=0.0, curvature_at_median=float(J.loc[(J['mod']-J['mod'].median()).abs().idxmin(), 'curvature']))

jn1 = jn_other(i_use, "ai_use")   # AI专项评估：在一般评估=0时，AI专项在何处使曲率显著<0
jn2 = jn_other(i_gen, "ai_gen")   # 一般网络评估：在AI专项=0时
print(f"  AI专项评估: 曲率显著<0 区域 {jn1['sig_neg_range']} (占网格 {jn1['frac_grid_sig']*100:.0f}%)", flush=True)
print(f"  一般网络评估: 曲率显著<0 区域 {jn2['sig_neg_range']} (占网格 {jn2['frac_grid_sig']*100:.0f}%)", flush=True)
pd.DataFrame([jn1, jn2]).to_csv(f"{OUT}/JohnsonNeyman_summary.csv", index=False, encoding="utf-8")

# ================================================================ (D) FDR 校正
print("\n[D] FDR (Benjamini–Hochberg) 多重比较校正", flush=True)
def bh(pvals):
    p = np.asarray(pvals, float); n = len(p); order = np.argsort(p)
    ranked = p[order]*n/np.arange(1, n+1)
    q_sorted = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n); q[order] = np.minimum(q_sorted, 1)
    return q

fdr_rows = []
task_csv = pd.read_csv("r_outputs_real/supplementary/ANA02_H1a_task_specific.csv")
b2 = task_csv[task_csv.term == "ai2"].copy()
b2["subset"] = b2["subset"].astype(str)
for p_val, label in zip(b2["p"].values, (b2["task"]+"_"+b2["subset"]).values):
    fdr_rows.append(dict(family="ANA02 four tasks x two samples", test=label, p=p_val))
meta_csv = pd.read_csv("r_outputs_real/meta_regression/univariate.csv")
for p_val, label in zip(meta_csv["p"].values, meta_csv["var"].values):
    fdr_rows.append(dict(family="meta-regression univariate", test=label, p=p_val))
FDR = pd.DataFrame(fdr_rows)
for fam, g in FDR.groupby("family"):
    q = bh(g["p"].values)
    FDR.loc[g.index, "q_BH"] = q
    FDR.loc[g.index, "sig_q05"] = q < 0.05
FDR.to_csv(f"{OUT}/FDR_BH.csv", index=False, encoding="utf-8")
print(FDR.to_string(index=False, float_format=lambda x: f"{x:.3g}"), flush=True)

# ================================================================ 汇总
summary = dict(
    N_pooled=len(d_main), N_H2b=len(dc),
    lm_verdicts=LM[["label", "verdict"]].to_dict("records"),
    best_forms={s: g.loc[g.AIC.idxmin(), "model"] for s, g in FORMS.groupby("sample", sort=False)},
    jn=[jn1, jn2],
    fdr_sig=int(FDR.sig_q05.sum()), fdr_total=len(FDR),
    sec=time.time()-t0)
with open(f"{OUT}/summary.json", "w", encoding="utf-8") as f:
    json.dump(summary, f, ensure_ascii=False, indent=1, default=float)
print(f"\nDONE sec={time.time()-t0:.1f}  OUT={OUT}", flush=True)
