# -*- coding: utf-8 -*-
"""
G04-01 取证：为"倒U假设与真实数据冲突"的处置提供客观证据（不替作者决定叙事）。
 证据A 分国家/经济体拟合：倒U(ai2<0显著)在多少子群成立？
 证据B 方向性混淆：AI 使用是否与学业困难(低ESCS)混淆？
 证据C 按 AI 使用频率分组的成绩与ESCS描述。
探索性取证：分国家用 5PV 均值、常规 WLS SE（非 BRR），仅供定性判断，非最终估计。
"""
import pyreadstat, numpy as np, pandas as pd
from scipy import stats
import os, time
t0=time.time()
OUT="r_outputs_real/g04_forensics"; os.makedirs(OUT,exist_ok=True)

cols=["CNT","ST438Q01DA","ST438Q02DA","ST438Q03DA","ST438Q04DA","AIUSESCH",
      "ESCS","ST004D01T","PV1SCIE","PV2SCIE","PV3SCIE","PV4SCIE","PV5SCIE","W_FSTUWT"]
print("reading...",flush=True)
df,_=pyreadstat.read_sav("data/CY09_MS_STU_PUF.sav",usecols=cols,apply_value_formats=False)
print(f"N={len(df):,}",flush=True)
for c in ["ST438Q01DA","ST438Q02DA","ST438Q03DA","ST438Q04DA"]:
    df.loc[df[c].isin([95,96,97,98,99]),c]=np.nan
df.loc[df["ESCS"].isin([95,96,97,98,99]),"ESCS"]=np.nan
df.loc[df["ST004D01T"].isin([8,9]),"ST004D01T"]=np.nan
df.loc[(~np.isfinite(df.AIUSESCH))|(df.AIUSESCH.abs()>10),"AIUSESCH"]=np.nan
pv=[f"PV{i}SCIE" for i in range(1,6)]
df["pvmean"]=df[pv].mean(axis=1)
df["gender"]=df["ST004D01T"]-1
df["ai4"]=df[["ST438Q01DA","ST438Q02DA","ST438Q03DA","ST438Q04DA"]].mean(axis=1,skipna=True)
d=df.dropna(subset=["AIUSESCH","ai4","ESCS","gender","pvmean","W_FSTUWT"]).copy()
d["ai"]=d.AIUSESCH-d.AIUSESCH.mean(); d["ai2"]=d.ai**2
d["escs_c"]=d.ESCS-d.ESCS.mean()
print(f"完整案例 N={len(d):,}",flush=True)

def wls(X,y,w):
    sw=np.sqrt(w); Xw=X*sw[:,None]; yw=y*sw
    b,*_=np.linalg.lstsq(Xw,yw,rcond=None)
    r=y-X@b
    meat=(X*w[:,None]*r[:,None]).T@(X*w[:,None]*r[:,None])
    inv=np.linalg.pinv(X.T@(X*w[:,None])); cov=inv@meat@inv
    se=np.sqrt(np.clip(np.diag(cov),0,None))
    return b,se

# ---- 证据A：分国家拟合 ------------------------------------------------------
rows=[]
for ctry,g in d.groupby("CNT"):
    n=len(g)
    if n<200: continue
    X=np.column_stack([np.ones(n),g.ai,g.ai2,g.escs_c,g.gender])
    w=g.W_FSTUWT.values/g.W_FSTUWT.mean()
    b,se=wls(X,g.pvmean.values,w)
    t=b[2]/se[2] if se[2]>0 else 0
    p=2*stats.t.sf(abs(t),max(n-5,1))
    rows.append(dict(country=ctry,n=n,ai2=b[2],se=se[2],p=p,
                     shape=("倒U" if (b[2]<0 and p<0.05) else ("U型" if (b[2]>0 and p<0.05) else "不显著")),
                     vertex=(-b[1]/(2*b[2]))+g.AIUSESCH.mean() if b[2]!=0 else np.nan))
A=pd.DataFrame(rows).sort_values("n",ascending=False)
A.to_csv(f"{OUT}/by_country.csv",index=False,encoding="utf-8")
print("\n===== 证据A：分国家倒U检验（n>=200）=====")
print(f"可检验国家数={len(A)}")
print(A["shape"].value_counts().to_string())
# p值合理性诊断：常规WLS的稳健SE在极端权重下可能失真，报告p<0.05占比
print(f"\np<0.05 的国家数: {(A.p<0.05).sum()} / {len(A)} （p值分布异常请谨慎解读）")
print("\n样本量最大的12国:")
print(A[["country","n","ai2","p","shape","vertex"]].head(12).to_string(index=False,float_format=lambda x:f"{x:.3f}"))
# 倒U国家中，顶点是否落在数据范围内
inv=A[A["shape"]=="倒U"].copy()
if len(inv):
    inr=[]
    for _,r in inv.iterrows():
        g=d[d.CNT==r.country]
        lo,hi=g.AIUSESCH.quantile(.01),g.AIUSESCH.quantile(.99)
        inr.append(lo<=r.vertex<=hi)
    inv["vertex_in_range"]=inr
    print(f"\n倒U国家中顶点落在数据范围内的: {int(inv.vertex_in_range.sum())}/{len(inv)}")

# ---- 证据B：方向性混淆（AI 使用 vs ESCS）-----------------------------------
print("\n===== 证据B：方向性混淆 =====")
r_all=np.corrcoef(d.AIUSESCH,d.ESCS)[0,1]
print(f"corr(AIUSESCH, ESCS) = {r_all:.3f}")
# 四题均值分5组
d["ai4grp"]=pd.qcut(d.ai4,5,labels=["Q1最低","Q2","Q3","Q4","Q5最高"])
B=d.groupby("ai4grp",observed=True).agg(n=("pvmean","size"),ai4=("ai4","mean"),
      escs=("ESCS","mean"),science=("pvmean","mean")).reset_index()
B.to_csv(f"{OUT}/by_ai_level.csv",index=False,encoding="utf-8")
print(B.to_string(index=False,float_format=lambda x:f"{x:.3f}"))

# ---- 证据C：整体二次项方向（引用主分析）------------------------------------
print("\n===== 证据C：整体 =====")
X=np.column_stack([np.ones(len(d)),d.ai,d.ai2,d.escs_c,d.gender])
w=norm=None
b,_=wls(X,d.pvmean.values,d.W_FSTUWT.values/d.W_FSTUWT.mean())
vtx=-b[1]/(2*b[2])+d.AIUSESCH.mean()
lo,hi=d.AIUSESCH.quantile(.01),d.AIUSESCH.quantile(.99)
print(f"整体 ai2={b[2]:.4f}（{'正/U型' if b[2]>0 else '负/倒U'}）；顶点原始尺度={vtx:.3f}；观测范围=[{lo:.3f},{hi:.3f}]；顶点在范围内={lo<=vtx<=hi}")
print(f"\n总耗时={time.time()-t0:.1f}s  输出={OUT}")
