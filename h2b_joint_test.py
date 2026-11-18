# -*- coding: utf-8 -*-
"""
H2b 同一模型系数差异检验：AI专项(ST436Q16DA) vs 一般网络(ST436Q02DA) 曲率交互
在同一回归中同时放入两个评估变量及其与 ai/ai2 的交互，
用 BRR 重复权重估计两个 ai2×eval 系数之差(contrast)的 SE 与 p。
判据(v2 H2b)：若 AI专项交互显著为负(强化倒U)且|系数|显著大于一般网络 → 支持AI特异性；
             若 AI专项不显著而一般网络显著 → 支持一般素养重解释。
方法与主分析一致：WLS + 逐PV + Rubin + Fay BRR(rho=0.5) + 国家内归一化权重。
"""
import pyreadstat, numpy as np, pandas as pd
from scipy import stats
import os, time
t0=time.time()
OUT="r_outputs_real/supplementary"; os.makedirs(OUT,exist_ok=True)
pv=[f"PV{i}SCIE" for i in range(1,6)]; rw=[f"W_FSTURWT{i}" for i in range(1,81)]
cols=["CNT","ST436Q16DA","ST436Q02DA","AIUSESCH","ESCS","ST004D01T"]+pv+["W_FSTUWT"]+rw
df,_=pyreadstat.read_sav("data/CY09_MS_STU_PUF.sav",usecols=cols,apply_value_formats=False)
for c in ["ST436Q16DA","ST436Q02DA"]: df.loc[df[c].isin([5,6,7,8,9]),c]=np.nan
df.loc[df["ESCS"].isin([95,96,97,98,99]),"ESCS"]=np.nan
df.loc[~df["ST004D01T"].isin([1,2]),"ST004D01T"]=np.nan
df.loc[(~np.isfinite(df["AIUSESCH"]))|(df["AIUSESCH"].abs()>10),"AIUSESCH"]=np.nan
for c in pv: df.loc[df[c]>=9997,c]=np.nan
df["gender"]=df["ST004D01T"]-1
d=df.dropna(subset=["AIUSESCH","ST436Q16DA","ST436Q02DA","ESCS","ST004D01T"]+pv+["W_FSTUWT"]+rw).copy()
d["ai"]=d["AIUSESCH"]-d["AIUSESCH"].mean(); d["ai2"]=d["ai"]**2
d["ai_use"]=d["ST436Q16DA"]-d["ST436Q16DA"].mean()
d["ai_gen"]=d["ST436Q02DA"]-d["ST436Q02DA"].mean()
d["escs_c"]=d["ESCS"]-d["ESCS"].mean()
X=pd.DataFrame({"const":1.0,"ai":d["ai"].values,"ai2":d["ai2"].values,
   "ai_use":d["ai_use"].values,"ai_gen":d["ai_gen"].values,
   "ai_x_use":d["ai"].values*d["ai_use"].values,"ai2_x_use":d["ai2"].values*d["ai_use"].values,
   "ai_x_gen":d["ai"].values*d["ai_gen"].values,"ai2_x_gen":d["ai2"].values*d["ai_gen"].values,
   "escs_c":d["escs_c"].values,"gender":d["gender"].values})
names=list(X.columns); Xm=X.values
i_use=names.index("ai2_x_use"); i_gen=names.index("ai2_x_gen")
def norm(w,g):
    s=pd.Series(w);gg=pd.Series(g);return w*s.groupby(gg).transform(lambda x:len(x)/x.sum()).values
w=norm(d["W_FSTUWT"].values,d["CNT"].values)
rwm=np.column_stack([norm(d[r].values,d["CNT"].values) for r in rw])
def wls(X,y,w):
    sw=np.sqrt(w);Xw=X*sw[:,None];yw=y*sw
    b,*_=np.linalg.lstsq(Xw,yw,rcond=None);return b
rho=0.5;R=len(rw)
# 逐PV：全样本系数 + BRR 重复系数（用于 contrast）
ests=[];contrasts=[];brr_all=[]
for pvc in pv:
    y=d[pvc].values.astype(float)
    b=wls(Xm,y,w);ests.append(b);contrasts.append(b[i_use]-b[i_gen])
    bs=np.empty((R,Xm.shape[1]))
    for r in range(R):
        wr=(1-rho)*w+rho*rwm[:,r]
        bs[r]=wls(Xm,y,wr)
    brr_all.append(bs)
# Rubin 合并主系数
def rubin(ests):
    Q=np.array(ests).T;M=Q.shape[1];qbar=Q.mean(1);b=((Q-qbar[:,None])**2).sum(1)/(M-1)
    return qbar,b
# BRR 方差：先在每个PV内算BRR协方差，再跨PV合并到 contrast
# contrast 的 BRR 变异 = mean over PV of Var_BRR(b_use - b_gen)
var_contrast_brr=[]
for bs in brr_all:
    dc=bs[:,i_use]-bs[:,i_gen]
    var_contrast_brr.append(((dc-dc.mean())**2).sum()/(R*(1-rho)**2))
U_c=np.mean(var_contrast_brr)
# 跨PV的contrast方差 b_c
cs=np.array(contrasts);cbar=cs.mean();B_c=cs.var(ddof=1)
T_c=U_c+(1+1/5)*B_c
se_c=np.sqrt(T_c);t_c=cbar/se_c;p_c=2*stats.t.sf(abs(t_c),df=1e6)
qbar,between=rubin(ests)
# 各系数BRR SE（首个PV的BRR，跨PV取均值方差）
var_b=np.mean([((bs-bs.mean(0))**2).sum(0)/(R*(1-rho)**2) for bs in brr_all],axis=0)
se_b=np.sqrt(var_b)
print(f"N={len(d):,}")
for nm in ["ai2","ai2_x_use","ai2_x_gen"]:
    i=names.index(nm);t=qbar[i]/se_b[i]
    print(f"{nm:12s} est={qbar[i]:+.4f} BRR_SE={se_b[i]:.4f} p={2*stats.t.sf(abs(t),1e6):.3g}")
print(f"CONTRAST (ai2_x_use - ai2_x_gen) = {cbar:+.4f}  BRR_SE={se_c:.4f}  p={p_c:.4g}")
rows=[dict(term=nm,est=qbar[i],brr_se=se_b[i],p=2*stats.t.sf(abs(qbar[i]/se_b[i]),1e6)) for i,nm in enumerate(names)]
rows.append(dict(term="CONTRAST_use_minus_gen",est=cbar,brr_se=se_c,p=p_c))
pd.DataFrame(rows).to_csv(f"{OUT}/ANA01_H2b_joint_contrast.csv",index=False,encoding="utf-8")
print(f"DONE sec={time.time()-t0:.1f}")
