# -*- coding: utf-8 -*-
"""
补充分析：补齐 v2 设计中三项此前完全未动的关键检验
  ANA-01 H2b 证伪对照：AI专项评估(ST436Q16DA) vs 一般网络信息评估(ST436Q02DA)
  ANA-02 H1a 任务版：ST438 四任务分别的曲线 + 用户内子样本（检验"非使用者抬高顶点"）
  ANA-03 H3/H4：ESCS 分层 + 性别×SES（机会梯度与成绩曲线分层）

方法：与 run_real_analysis.py 完全一致——纯 NumPy WLS + 逐 PV 拟合 + Rubin 合并
      + Fay BRR(rho=0.5, 80 重复权重) SE + 国家内归一化权重。
立场：如实报告，不预设、不迎合论文叙事。未严格执行 v2 计划处均在报告注明。
"""
import pyreadstat, numpy as np, pandas as pd
from scipy import stats
import os, time, json
t0 = time.time()
OUT = "r_outputs_real/supplementary"
os.makedirs(OUT, exist_ok=True)

pv = [f"PV{i}SCIE" for i in range(1,6)]
rw = [f"W_FSTURWT{i}" for i in range(1,81)]
cols = ["CNT","ST436Q16DA","ST436Q02DA",
        "ST438Q01DA","ST438Q02DA","ST438Q03DA","ST438Q04DA",
        "ESCS","ST004D01T"] + pv + ["W_FSTUWT"] + rw
print("reading sav...", flush=True)
df, _ = pyreadstat.read_sav("data/CY09_MS_STU_PUF.sav", usecols=cols, apply_value_formats=False)
print(f"N={len(df):,} read_sec={time.time()-t0:.1f}", flush=True)

# ---- 缺失重编码（严格按 codebook 值标签）----
for c in ["ST436Q16DA","ST436Q02DA"]:
    df.loc[df[c].isin([5,6,7,8,9]), c] = np.nan
for c in ["ST438Q01DA","ST438Q02DA","ST438Q03DA","ST438Q04DA"]:
    df.loc[df[c].isin([95,96,97,98,99]), c] = np.nan
df.loc[df["ESCS"].isin([95,96,97,98,99]), "ESCS"] = np.nan
df.loc[~df["ST004D01T"].isin([1,2]), "ST004D01T"] = np.nan
for c in pv:
    df.loc[df[c]>=9997, c] = np.nan
df["gender"] = df["ST004D01T"] - 1     # 0=female, 1=male

def norm_weights(w, countries):
    s = pd.Series(w); g = pd.Series(countries)
    scale = s.groupby(g).transform(lambda x: len(x)/x.sum()).values
    return w*scale

def wls(X, y, w):
    sw = np.sqrt(w); Xw = X*sw[:,None]; yw = y*sw
    beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    resid = y - X@beta
    Z = X*w[:,None]*resid[:,None]
    meat = Z.T@Z
    inv = np.linalg.pinv(X.T@(X*w[:,None]))
    cov = inv@meat@inv
    return beta, np.sqrt(np.clip(np.diag(cov),0,None))

def rubin(ests, ses):
    Q=np.array(ests).T; U=np.array(ses).T**2; M=Q.shape[1]
    qbar=Q.mean(1); ubar=U.mean(1)
    b=((Q-qbar[:,None])**2).sum(1)/(M-1)
    T=ubar+(1+1/M)*b
    dfree=(M-1)*(1+ubar/((1+1/M)*np.where(b>0,b,np.nan)))**2
    dfree[~np.isfinite(dfree)|(b<=0)]=1e6
    t=qbar/np.sqrt(np.where(T>0,T,np.nan)); p=2*stats.t.sf(np.abs(t),dfree)
    return qbar,np.sqrt(T),t,p,dfree

rho=0.5
def brr_se(X, y, w_full, rw_mat):
    R=rw_mat.shape[1]
    bs=np.empty((R,X.shape[1]))
    for r in range(R):
        wr=(1-rho)*w_full+rho*rw_mat[:,r]
        bs[r],_=wls(X,y,wr)
    return np.sqrt(((bs-bs.mean(0))**2).sum(0)/(R*(1-rho)**2))

def fit_family(X, names, dsub, label):
    """对给定设计矩阵，逐 PV 拟合 + BRR SE + Rubin 合并；返回系数表 DataFrame"""
    w = norm_weights(dsub["W_FSTUWT"].values, dsub["CNT"].values)
    rwm = np.column_stack([norm_weights(dsub[r].values, dsub["CNT"].values) for r in rw])
    ests=[]; ses=[]
    for pvc in pv:
        y=dsub[pvc].values.astype(float)
        b,_=wls(X,y,w)
        ests.append(b); ses.append(brr_se(X,y,w,rwm))
    q,se,t,p,_=rubin(ests,ses)
    tab=pd.DataFrame(dict(term=names,est=q,se=se,t=t,p=p))
    tab["model"]=label
    return tab

def make_sub(cols_need):
    d=df.dropna(subset=cols_need+pv+["W_FSTUWT"]+rw).copy()
    d["escs_c"]=d["ESCS"]-d["ESCS"].mean()
    return d

results={}

# ============ ANA-01: H2b 证伪对照 ============
print("\n[ANA-01] H2b AI专项 vs 一般网络评估", flush=True)
need01=["ST436Q16DA","ST436Q02DA","ESCS","ST004D01T","AIUSESCH"] if "AIUSESCH" in df else ["ST436Q16DA","ST436Q02DA","ESCS","ST004D01T"]
# 主分析用 AIUSESCH 作为 ai 测量；本脚本未读 AIUSESCH，补读其所需——用四题均值近似会偏离主分析。
# 为保持一致，重新从已存 analysis_ready.csv 不可行（缺新变量），故补读 AIUSESCH。
print("  补读 AIUSESCH 以复用主分析测量...", flush=True)
extra,_=pyreadstat.read_sav("data/CY09_MS_STU_PUF.sav",usecols=["CNTSTUID","AIUSESCH"],apply_value_formats=False)
df["AIUSESCH"]=extra["AIUSESCH"].values
df.loc[(~np.isfinite(df["AIUSESCH"]))|(df["AIUSESCH"].abs()>10),"AIUSESCH"]=np.nan

d01=make_sub(["AIUSESCH","ST436Q16DA","ST436Q02DA","ESCS","ST004D01T"])
d01["ai"]=d01["AIUSESCH"]-d01["AIUSESCH"].mean()
d01["ai2"]=d01["ai"]**2
d01["ai_use"]=d01["ST436Q16DA"]-d01["ST436Q16DA"].mean()   # AI专项评估
d01["ai_gen"]=d01["ST436Q02DA"]-d01["ST436Q02DA"].mean()   # 一般网络评估
d01["gen_c"]=d01["ESCS"]-d01["ESCS"].mean()
tabs01=[]
for tag,var in [("AI_specific(ST436Q16DA)","ai_use"),("General_online(ST436Q02DA)","ai_gen")]:
    X=pd.DataFrame({"const":1.0,"ai":d01["ai"].values,"ai2":d01["ai2"].values,
                    var:d01[var].values,
                    f"ai_x_{var}":d01["ai"].values*d01[var].values,
                    f"ai2_x_{var}":d01["ai2"].values*d01[var].values,
                    "escs_c":d01["gen_c"].values,"gender":d01["gender"].values})
    t=fit_family(X.values,list(X.columns),d01,f"H2b_{tag}")
    tabs01.append(t)
    print(f"  {tag}: N={len(d01):,}", flush=True)
    print(t[["term","est","se","p"]].to_string(index=False,float_format=lambda x:f"{x:.4g}"),flush=True)
res01=pd.concat(tabs01,ignore_index=True)
res01.to_csv(f"{OUT}/ANA01_H2b_falsification.csv",index=False,encoding="utf-8")
results["ANA01_N"]=len(d01)

# ============ ANA-02: H1a 任务版（全样本 + 用户内） ============
print("\n[ANA-02] H1a 四任务曲线（全样本 vs 用户内）", flush=True)
task_cols=["ST438Q01DA","ST438Q02DA","ST438Q03DA","ST438Q04DA"]
task_names=["summarise","prelim_research","drafting","help_learn"]
tabs02=[]
for tc,tn in zip(task_cols,task_names):
    for subset in ["all","users_only"]:
        d02=make_sub([tc,"ESCS","ST004D01T"])
        if subset=="users_only":
            d02=d02[d02[tc]>1].copy()   # 删除"从不/几乎从不"(=1)
        d02["ai"]=d02[tc]-d02[tc].mean()
        d02["ai2"]=d02["ai"]**2
        d02["escs_c"]=d02["ESCS"]-d02["ESCS"].mean()
        X=pd.DataFrame({"const":1.0,"ai":d02["ai"].values,"ai2":d02["ai2"].values,
                        "escs_c":d02["escs_c"].values,"gender":d02["gender"].values})
        t=fit_family(X.values,list(X.columns),d02,f"{tn}_{subset}")
        b1=float(t.loc[t.term=="ai","est"].iloc[0]); b2=float(t.loc[t.term=="ai2","est"].iloc[0])
        p2=float(t.loc[t.term=="ai2","p"].iloc[0])
        vertex=-b1/(2*b2) if b2!=0 else np.nan
        # 顶点是否落在该任务原始1-5范围内（中心化还原：vertex_raw = vertex + mean）
        vraw=vertex+d02[tc].mean() if np.isfinite(vertex) else np.nan
        t["vertex_centered"]=vertex; t["vertex_raw"]=vraw
        t["vertex_in_1_5"]=(vraw>=1)&(vraw<=5)
        t["N"]=len(d02); t["task"]=tn; t["subset"]=subset
        tabs02.append(t)
        print(f"  {tn:16s} {subset:10s} N={len(d02):,} ai2={b2:+.3f} p={p2:.4g} vertex_raw={vraw:.2f} in_range={t['vertex_in_1_5'].iloc[0]}",flush=True)
res02=pd.concat(tabs02,ignore_index=True)
res02.to_csv(f"{OUT}/ANA02_H1a_task_specific.csv",index=False,encoding="utf-8")

# ============ ANA-03: H3/H4 ESCS分层 + 性别×SES ============
print("\n[ANA-03] H3/H4 ESCS分层 + 性别×SES", flush=True)
# H3 机会梯度：以 AI 专项评估机会(ST436Q16DA)为因变量，WLS 近似（v2 计划为有序logit，见报告注明）
d03=make_sub(["ST436Q16DA","ST436Q02DA","ESCS","ST004D01T","AIUSESCH"])
d03["ai"]=d03["AIUSESCH"]-d03["AIUSESCH"].mean(); d03["ai2"]=d03["ai"]**2
d03["escs_c"]=d03["ESCS"]-d03["ESCS"].mean()
# 3a 机会梯度（因变量=评估机会）
tabs03=[]
for yv,ylab in [("ST436Q16DA","AI_eval_opportunity"),("ST436Q02DA","General_eval_opportunity")]:
    X=pd.DataFrame({"const":1.0,"escs_c":d03["escs_c"].values,"gender":d03["gender"].values,
                    "escs_x_gender":d03["escs_c"].values*d03["gender"].values,
                    "ai":d03["ai"].values})
tabs03=[]
for yv,ylab in [("ST436Q16DA","AI_eval_opportunity"),("ST436Q02DA","General_eval_opportunity")]:
    Xd=pd.DataFrame({"const":1.0,"escs_c":d03["escs_c"].values,"gender":d03["gender"].values,
                     "escs_x_gender":d03["escs_c"].values*d03["gender"].values,"ai":d03["ai"].values})
    w=norm_weights(d03["W_FSTUWT"].values,d03["CNT"].values)
    rwm=np.column_stack([norm_weights(d03[r].values,d03["CNT"].values) for r in rw])
    yv_mean=d03[yv].values.astype(float)
    # 因变量不是 PV，用单次 WLS + BRR SE（pv 列表不适用），单独处理
    b,_=wls(Xd.values,yv_mean,w)
    R=rwm.shape[1]; bs=np.empty((R,Xd.shape[1]))
    for r in range(R):
        wr=(1-rho)*w+rho*rwm[:,r]
        bs[r],_=wls(Xd.values,yv_mean,wr)
    se=np.sqrt(((bs-bs.mean(0))**2).sum(0)/(R*(1-rho)**2))
    tstat=b/np.where(se>0,se,np.nan); pp=2*stats.t.sf(np.abs(tstat),df=max(len(d03)-Xd.shape[1],1))
    t=pd.DataFrame(dict(term=list(Xd.columns),est=b,se=se,t=tstat,p=pp))
    t["model"]=f"H3_{ylab}"; t["N"]=len(d03)
    tabs03.append(t)
    print(f"  H3 {ylab}: N={len(d03):,} escs_c b={b[1]:+.4f} p={pp[1]:.4g}",flush=True)

# 3b 成绩曲线按 ESCS 三分层（国家内分位）
d03["escs_tier"]=d03.groupby("CNT")["ESCS"].transform(lambda x: pd.qcut(x,3,labels=[0,1,2]) if x.notna().sum()>30 else np.nan)
tabs3b=[]
for tier,lab in [(0,"low_ESCS"),(1,"mid_ESCS"),(2,"high_ESCS")]:
    ds=d03[d03["escs_tier"]==tier].copy()
    if len(ds)<1000: continue
    ds["ai"]=ds["AIUSESCH"]-ds["AIUSESCH"].mean(); ds["ai2"]=ds["ai"]**2
    ds["escs_c"]=ds["ESCS"]-ds["ESCS"].mean()
    X=pd.DataFrame({"const":1.0,"ai":ds["ai"].values,"ai2":ds["ai2"].values,
                    "escs_c":ds["escs_c"].values,"gender":ds["gender"].values})
    t=fit_family(X.values,list(X.columns),ds,f"H4_{lab}")
    b1=float(t.loc[t.term=="ai","est"].iloc[0]);b2=float(t.loc[t.term=="ai2","est"].iloc[0])
    t["vertex_centered"]=-b1/(2*b2) if b2!=0 else np.nan
    t["N"]=len(ds)
    tabs3b.append(t)
    print(f"  H4 {lab}: N={len(ds):,} ai={b1:+.3f} ai2={b2:+.3f} p={float(t.loc[t.term=='ai2','p'].iloc[0]):.4g}",flush=True)
res03=pd.concat(tabs03+tabs3b,ignore_index=True)
res03.to_csv(f"{OUT}/ANA03_H3H4_stratification.csv",index=False,encoding="utf-8")

summary=dict(ANA01_N=len(d01),ANA02_tasks=task_names,ANA03_N=len(d03),
             total_sec=round(time.time()-t0,1))
with open(f"{OUT}/summary.json","w",encoding="utf-8") as f: json.dump(summary,f,ensure_ascii=False,indent=1)
print(f"\nDONE total_sec={time.time()-t0:.1f} OUT={OUT}",flush=True)
