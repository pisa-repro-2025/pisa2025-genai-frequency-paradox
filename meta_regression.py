# -*- coding: utf-8 -*-
"""
S-02：跨国异质性解释 —— 元回归(meta-regression) + Egger 检验
目标：解释 formal_by_country.py 中 ai2 系数的 I²=97% 异质性来源。
因变量：各国 ai2（倒U曲率）; 自变量：国家层面特征。
方法：随机效应元回归（逆方差权重 1/(se²+tau²)，DL 型残差异质性迭代2轮）
      + 单变量筛选 + 多变量模型 + Egger 回归检验小样本偏倚。
立场：如实报告。R² 若很低则明说"国家层面特征解释力有限"，不强行叙事。
"""
import pyreadstat, numpy as np, pandas as pd
from scipy import stats
import os, time, warnings
warnings.filterwarnings("ignore")
t0 = time.time()
OUT = "r_outputs_real/meta_regression"; os.makedirs(OUT, exist_ok=True)

# 各国 ai2 与 se
A = pd.read_csv("r_outputs_real/formal_by_country/by_country_formal.csv")
A = A[A["ai2"].notna() & A["ai2_se"].notna() & (A["ai2_se"]>0)].copy().reset_index(drop=True)

# 从 sav 聚合国家层面协变量（用归一化主权重）
need = ["CNT","ESCS","AIUSESCH","ST437Q01DA","ST437Q02DA",
        "ST004D01T","W_FSTUWT"] + [f"PV{i}SCIE" for i in range(1,6)]
print("reading sav for country covariates...", flush=True)
df, _ = pyreadstat.read_sav("data/CY09_MS_STU_PUF.sav", usecols=need, apply_value_formats=False)
pv = [f"PV{i}SCIE" for i in range(1,6)]
for c in ["ST437Q01DA","ST437Q02DA"]:
    df.loc[df[c].isin([5,6,7,8,9]), c] = np.nan
df.loc[df["ESCS"].isin([95,96,97,98,99]), "ESCS"] = np.nan
df.loc[df["ST004D01T"].isin([8,9]), "ST004D01T"] = np.nan
df.loc[(~np.isfinite(df.AIUSESCH))|(df.AIUSESCH.abs()>10), "AIUSESCH"] = np.nan

def wavg(x, w):
    m = np.isfinite(x) & np.isfinite(w) & (w>0)
    return np.sum(x[m]*w[m])/np.sum(w[m]) if m.sum() else np.nan

rows=[]
for ctry, g in df.groupby("CNT"):
    g = g[g["W_FSTUWT"]>0]
    if len(g) < 200: continue
    w = g["W_FSTUWT"].values/g["W_FSTUWT"].mean()
    aiopp = np.nanmean(np.column_stack([g["ST437Q01DA"].values,g["ST437Q02DA"].values]),axis=1)
    rows.append(dict(country=ctry,
        cnt_escs=wavg(g["ESCS"].values,w),
        cnt_aiuse=wavg(g["AIUSESCH"].values,w),
        cnt_science=wavg(g[pv].mean(axis=1).values,w),
        cnt_aiopp=wavg(aiopp,w),
        cnt_female=wavg((g["ST004D01T"]==2).astype(float).values,w)))
C = pd.DataFrame(rows)
M = A.merge(C, on="country", how="inner").dropna()
print("merge 后国家数 k=%d" % len(M), flush=True)

def dl_tau2(y, v):
    k=len(y); w=1/v; mu=(w*y).sum()/w.sum()
    Q=(w*(y-mu)**2).sum(); Csum=w.sum()-(w**2).sum()/w.sum()
    return max(0.0,(Q-(k-1))/Csum) if Csum>0 else 0.0, Q

TAU2_TOTAL, Q_TOTAL = dl_tau2(M.ai2.values, M.ai2_se.values**2)

def meta_reg(Xnames, dfm):
    yv = dfm["ai2"].values.astype(float)
    v  = dfm["ai2_se"].values.astype(float)**2
    Xm = np.column_stack([np.ones(len(dfm))]+[dfm[n].values.astype(float) for n in Xnames])
    k = len(yv); p = Xm.shape[1]
    tau2 = TAU2_TOTAL                      # 锁定为无条件 DL 估计（标准做法）
    w = 1/(v+tau2); sw = np.sqrt(w)
    mu_re = (w*yv).sum()/w.sum()
    Qe_null = (w*(yv-mu_re)**2).sum()      # 空模型(仅截距)在RE权重下的异质性
    beta,*_ = np.linalg.lstsq(Xm*sw[:,None], yv*sw, rcond=None)
    cov = np.linalg.pinv((Xm*w[:,None]).T@Xm)
    se = np.sqrt(np.clip(np.diag(cov),0,None))
    t = beta/np.where(se>0,se,np.nan); pv = 2*stats.t.sf(np.abs(t),max(k-p,1))
    resid = yv - Xm@beta
    Qe = (w*resid**2).sum()                # 残差异质性 Q_E ~ chi2(k-p)
    I2 = max(0.0,(Qe-(k-p))/Qe)*100 if Qe>0 else 0.0
    r2 = max(0.0,(Qe_null-Qe)/Qe_null)*100 if Qe_null>0 else 0.0  # 同一RE权重下解释占比
    return beta,se,pv,tau2,I2,k,Qe,r2

cands=["cnt_escs","cnt_science","cnt_aiuse","cnt_aiopp","cnt_female"]
single=[]
for c in cands:
    dsub=M[["ai2","ai2_se",c]].dropna()
    if len(dsub)<30: continue
    beta,se,pv,tau2,I2,k,Qe,r2=meta_reg([c],dsub)
    single.append(dict(var=c,k=k,beta=beta[1],se=se[1],p=pv[1],
                       I2_resid=I2,Qe=Qe,R2_explained=r2))
S=pd.DataFrame(single).sort_values("p")
S.to_csv(f"{OUT}/univariate.csv",index=False,encoding="utf-8")
print("\n基线: tau2_DL=%.3f  Q_total=%.1f  I2_total=%.1f%%"%(TAU2_TOTAL,Q_TOTAL,max(0,(Q_TOTAL-(len(M)-1))/Q_TOTAL)*100),flush=True)
print("\n===== 单变量元回归(按p排序) =====",flush=True)
print(S.to_string(index=False,float_format=lambda x:"%.4g"%x),flush=True)

# 多变量模型：样本量允许下取最多3个变量（避免过拟合 k/p>=10）
sel=S["var"].head(3).tolist() if len(S)>=3 else cands[:2]
dsub=M[["ai2","ai2_se"]+sel].dropna()
beta,se,pv,tau2,I2,k,Qe,r2=meta_reg(sel,dsub)
names=["const"]+sel
mult=pd.DataFrame(dict(term=names,beta=beta,se=se,p=pv))
mult.to_csv(f"{OUT}/multivariate.csv",index=False,encoding="utf-8")
print("\n===== 多变量元回归: %s ====="%sel,flush=True)
print("k=%d tau2=%.3f I2_resid=%.1f%% Qe=%.1f R2_explained=%.1f%%"%(k,tau2,I2,Qe,r2),flush=True)
print(mult.to_string(index=False,float_format=lambda x:"%.4g"%x),flush=True)

# Egger 回归：ai2/se ~ 1/se (截距检验不对称)
yy=M.ai2.values/M.ai2_se.values
xx=1/M.ai2_se.values
sl,ic,r,pv_e,se_e=stats.linregress(xx,yy)
print("\n===== Egger 回归(小样本/发表偏倚) =====",flush=True)
print("截距=%.4f p=%.4g (截距≠0 提示不对称)"%(ic,pv_e),flush=True)

# 汇总
with open(f"{OUT}/summary.txt","w",encoding="utf-8") as f:
    f.write("S-02 元回归: k=%d 国家, 因变量=各国ai2\n"%len(M))
    f.write("基线 tau2_DL=%.3f Q_total=%.1f I2_total=%.1f%%\n"%(TAU2_TOTAL,Q_TOTAL,max(0,(Q_TOTAL-(len(M)-1))/Q_TOTAL)*100))
    f.write("单变量显著(p<0.05): %s\n"%(S[S.p<0.05]["var"].tolist() if len(S) else "无"))
    f.write("多变量模型: %s; 残差I2=%.1f%%; R2_explained=%.1f%%\n"%(sel,I2,r2))
    f.write("Egger 截距=%.4f p=%.4g\n"%(ic,pv_e))
print("\n总耗时=%.1fs 输出=%s"%(time.time()-t0,OUT),flush=True)
