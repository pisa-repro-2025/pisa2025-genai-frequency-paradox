# -*- coding: utf-8 -*-
"""
G04-01 正式估计：分国家/经济体 倒U型 逐PV + BRR 重复权重 + Rubin合并 + 跨国元分析
用途：把上一轮 G04_forensics.py 的定性取证（PV均值+常规SE）升级为可写入论文
      "结果"章节的正式估计。默认叙事方向=分国家异质倒U（数据支持最强）。
立场：如实报告；不预设每个国家都倒U，I² 异质性本身即是结论。

方法：
  对每个国家(n>=200)：
    对每个 PV(1-5)：
      主权重(国家内归一化) WLS: pv ~ 1 + ai + ai2 + escs_c + gender
      Fay BRR(rho=0.5, 80个重复权重, 国家内归一化) 估计 SE
    5个PV用 Rubin 规则合并 -> ai(线性)/ai2(二次) 的 est, se, p, 顶点
    判定 shape: 倒U(ai2<0,p<.05) / U型(ai2>0,p<.05) / 不显著
    顶点落在该国 AIUSESCH 1-99% 分位内 -> vertex_in_range
  跨国元分析(基于各国 ai2 的 est 与 se)：
    固定效应(逆方差) + 随机效应(DerSimonian-Laird) + Q/I^2/tau^2
"""
import pyreadstat, numpy as np, pandas as pd
from scipy import stats
import os, time, warnings
warnings.filterwarnings("ignore")
t0 = time.time()
OUT = "r_outputs_real/formal_by_country"
os.makedirs(OUT, exist_ok=True)

cols = ["CNT","AIUSESCH","ESCS","ST004D01T",
        "PV1SCIE","PV2SCIE","PV3SCIE","PV4SCIE","PV5SCIE",
        "W_FSTUWT"] + [f"W_FSTURWT{i}" for i in range(1,81)]
print("reading sav (usecols=%d)..." % len(cols), flush=True)
df, _ = pyreadstat.read_sav("data/CY09_MS_STU_PUF.sav", usecols=cols, apply_value_formats=False)
print("N=%d read_sec=%.1f" % (len(df), time.time()-t0), flush=True)

# 缺失值处理（与主分析一致）
df.loc[df["ESCS"].isin([95,96,97,98,99]), "ESCS"] = np.nan
df.loc[df["ST004D01T"].isin([8,9]), "ST004D01T"] = np.nan
df.loc[(~np.isfinite(df["AIUSESCH"])) | (df["AIUSESCH"].abs()>10), "AIUSESCH"] = np.nan
pv = [f"PV{i}SCIE" for i in range(1,6)]
rw = [f"W_FSTURWT{i}" for i in range(1,81)]
need = ["AIUSESCH","ESCS","ST004D01T","W_FSTUWT"] + rw + pv
d = df.dropna(subset=need).copy()
gm = d["AIUSESCH"].mean()          # 全局中心化常数（保证各国 ai2 可比）
d["gender"] = d["ST004D01T"] - 1
d["ai"]  = d["AIUSESCH"] - gm
d["ai2"] = d["ai"]**2
d["escs_c"] = d["ESCS"] - d["ESCS"].mean()
print("完整案例 N=%d (%.1f%%) ai2与中心化常数无关, gm=%.3f" %
      (len(d), 100*len(d)/len(df), gm), flush=True)

# 国家内归一化权重（主权重与80个重复权重）
def norm_in_country(w, ctry):
    s = pd.Series(w); g = pd.Series(ctry)
    scale = s.groupby(g).transform(lambda x: len(x)/x.sum()).values
    return w*scale
d["w0"] = norm_in_country(d["W_FSTUWT"].values, d["CNT"].values)
for r in rw:
    d[r+"_n"] = norm_in_country(d[r].values, d["CNT"].values)
rw_n = [r+"_n" for r in rw]

def wls_beta(X, y, w):
    sw = np.sqrt(w); Xw = X*sw[:,None]; yw = y*sw
    b, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    return b

RHO = 0.5; R = 80
def brr_se(X, y, w0, rwmat):
    b_full = wls_beta(X, y, w0)
    bs = np.empty((R, X.shape[1]))
    for r in range(R):
        wr = (1-RHO)*w0 + RHO*rwmat[:,r]
        wr = wr/wr.mean()
        bs[r] = wls_beta(X, y, wr)
    var = ((bs-bs.mean(0))**2).sum(0)/(R*(1-RHO)**2)
    return b_full, np.sqrt(var)

def rubin(ests, ses):
    Q = np.asarray(ests).T; U = np.asarray(ses).T**2; M = Q.shape[1]
    qbar = Q.mean(1); ubar = U.mean(1)
    b = ((Q-qbar[:,None])**2).sum(1)/(M-1)
    T = ubar + (1+1.0/M)*b
    dfree = (M-1)*(1+ubar/((1+1.0/M)*np.where(b>0,b,np.nan)))**2
    dfree = np.where(np.isfinite(dfree) & (b>0), dfree, 1e6)
    t = qbar/np.sqrt(np.where(T>0,T,np.nan)); p = 2*stats.t.sf(np.abs(t), dfree)
    return qbar, np.sqrt(np.where(T>0,T,0)), p

rows = []
countries = sorted(d["CNT"].unique())
print("可检验国家(n>=200)遍历开始, 总计%d个国家/经济体..." % len(countries), flush=True)
done = 0
for ctry in countries:
    g = d[d["CNT"]==ctry]
    n = len(g)
    if n < 200:
        continue
    X = np.column_stack([np.ones(n), g["ai"].values, g["ai2"].values,
                         g["escs_c"].values, g["gender"].values])
    w0 = g["w0"].values
    rwmat = g[rw_n].values
    ests, ses = [], []
    for pvc in pv:
        yv = g[pvc].astype(float).values
        b, se = brr_se(X, yv, w0, rwmat)
        ests.append(b); ses.append(se)
    qbar, se, p = rubin(ests, ses)
    b1, b2 = qbar[1], qbar[2]
    se2 = se[2]; p2 = p[2]
    shape = "倒U" if (b2<0 and p2<0.05) else ("U型" if (b2>0 and p2<0.05) else "不显著")
    vtx = (-b1/(2*b2)) + gm if b2 != 0 else np.nan
    lo, hi = g["AIUSESCH"].quantile(.01), g["AIUSESCH"].quantile(.99)
    in_range = bool(lo <= vtx <= hi)
    rows.append(dict(country=ctry, n=n, ai=b1, ai_se=se[1], ai_p=p[1],
                     ai2=b2, ai2_se=se2, ai2_p=p2, shape=shape,
                     vertex=vtx, vlo=lo, vhi=hi, vertex_in_range=in_range))
    done += 1
    if done % 10 == 0:
        print("  已处理 %d 国, 用时 %.1fs" % (done, time.time()-t0), flush=True)

A = pd.DataFrame(rows).sort_values("n", ascending=False)
A.to_csv(f"{OUT}/by_country_formal.csv", index=False, encoding="utf-8")
print("\n===== 正式分国估计: shape 分布 =====")
print(A["shape"].value_counts().to_string(), flush=True)
inv = A[A["shape"]=="倒U"]
if len(inv):
    print("倒U国家中顶点落在数据范围内: %d/%d" %
          (int(inv["vertex_in_range"].sum()), len(inv)), flush=True)

# ---- 跨国元分析(基于各国 ai2) ----------------------------------------------
print("\n===== 跨国元分析(ai2 二次项系数) =====", flush=True)
b = A["ai2"].values.astype(float); se = A["ai2_se"].values.astype(float)
ok = np.isfinite(b) & np.isfinite(se) & (se>0)
b, se, ctry_v = b[ok], se[ok], A["country"].values[ok]
k = len(b)
w_fix = 1/se**2
b_fe = (w_fix*b).sum()/w_fix.sum(); se_fe = 1/np.sqrt(w_fix.sum())
Q = (w_fix*(b-b_fe)**2).sum()
C = w_fix.sum() - (w_fix**2).sum()/w_fix.sum()
tau2 = max(0.0, (Q-(k-1))/C)
I2 = max(0.0, (Q-(k-1))/Q)*100 if Q>0 else 0.0
w_re = 1/(se**2+tau2)
b_re = (w_re*b).sum()/w_re.sum(); se_re = 1/np.sqrt(w_re.sum())
p_fe = 2*stats.norm.sf(abs(b_fe/se_fe)); p_re = 2*stats.norm.sf(abs(b_re/se_re))
meta = pd.DataFrame([
    dict(model="固定效应", k=k, est=b_fe, se=se_fe, p=p_fe, tau2=np.nan, I2=I2, Q=Q),
    dict(model="随机效应DL", k=k, est=b_re, se=se_re, p=p_re, tau2=tau2, I2=I2, Q=Q),
])
meta.to_csv(f"{OUT}/meta_analysis.csv", index=False, encoding="utf-8")
print(meta.to_string(index=False, float_format=lambda x:"%.4f"%x), flush=True)
print("\n解读: ai2<0=倒U, ai2>0=U型; I2高=跨国高度异质(pooled单一结论不成立)", flush=True)
print("总耗时=%.1fs 输出=%s" % (time.time()-t0, OUT), flush=True)
