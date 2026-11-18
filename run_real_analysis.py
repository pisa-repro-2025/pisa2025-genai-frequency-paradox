# -*- coding: utf-8 -*-
"""
真实 PISA 2025 数据分析：AI 使用频率 x 科学成绩倒U型 + AI素养教学调节
决策依据 DECISIONS.md。方法：纯 NumPy 加权最小二乘(WLS) + 逐PV拟合 + Rubin合并
       + Fay BRR 重复权重估计抽样方差。
立场：如实报告，倒U是否成立以数据为准，不预设、不迎合论文叙事。
说明：本机安全策略拦截所有含DLL的R包(statsmodels/survey/data.table等)，
      故用 numpy 手写 WLS（数学等价），不依赖被拦截的编译扩展。
"""
import pyreadstat, numpy as np, pandas as pd
from scipy import stats
import os, time
t0 = time.time()
OUT = "r_outputs_real/real_results"
os.makedirs(OUT, exist_ok=True)

cols = ["CNT","CNTSCHID",
        "ST438Q01DA","ST438Q02DA","ST438Q03DA","ST438Q04DA","AIUSESCH",
        "ST436Q02DA","ESCS","ST004D01T",
        "PV1SCIE","PV2SCIE","PV3SCIE","PV4SCIE","PV5SCIE",
        "W_FSTUWT"] + [f"W_FSTURWT{i}" for i in range(1,81)]
print("reading sav (usecols=96)...", flush=True)
df, _ = pyreadstat.read_sav("data/CY09_MS_STU_PUF.sav", usecols=cols, apply_value_formats=False)
print(f"N={len(df):,}  read_sec={time.time()-t0:.1f}", flush=True)

# 缺失值重编码
for c in ["ST438Q01DA","ST438Q02DA","ST438Q03DA","ST438Q04DA"]:
    df.loc[df[c].isin([95,96,97,98,99]), c] = np.nan
df.loc[df["ST436Q02DA"].isin([5,6,7,8,9]), "ST436Q02DA"] = np.nan
df.loc[df["ESCS"].isin([95,96,97,98,99]), "ESCS"] = np.nan
df.loc[df["ST004D01T"].isin([8,9]), "ST004D01T"] = np.nan
df.loc[(~np.isfinite(df["AIUSESCH"])) | (df["AIUSESCH"].abs()>10), "AIUSESCH"] = np.nan

df["gender"] = df["ST004D01T"] - 1
df["ai"]  = df["AIUSESCH"] - df["AIUSESCH"].mean()
df["ai2"] = df["ai"]**2
df["instr"] = df["ST436Q02DA"] - df["ST436Q02DA"].mean()
df["escs_c"] = df["ESCS"] - df["ESCS"].mean()

pv = [f"PV{i}SCIE" for i in range(1,6)]
rw = [f"W_FSTURWT{i}" for i in range(1,81)]
need = ["ai","ai2","instr","escs_c","gender","W_FSTUWT"] + rw
d = df.dropna(subset=need + pv).copy()
print(f"完整案例 N={len(d):,} ({100*len(d)/len(df):.1f}%)", flush=True)

# 国家内归一化权重
def norm_weights(w, countries):
    s = pd.Series(w); g = pd.Series(countries)
    scale = s.groupby(g).transform(lambda x: len(x)/x.sum()).values
    return w*scale
w_main = norm_weights(d["W_FSTUWT"].values, d["CNT"].values)
rw_norm = np.column_stack([norm_weights(d[r].values, d["CNT"].values) for r in rw])

def design(kind):
    base = pd.DataFrame({
        "const":1.0, "ai":d["ai"].values, "ai2":d["ai2"].values,
        "escs_c":d["escs_c"].values, "gender":d["gender"].values})
    if kind=="mod":
        base["instr"]=d["instr"].values
        base["ai_x_instr"]=d["ai"].values*d["instr"].values
        base["ai2_x_instr"]=d["ai2"].values*d["instr"].values
    return base.values, list(base.columns)

def wls(X, y, w):
    sw = np.sqrt(w)
    Xw = X*sw[:,None]; yw = y*sw
    beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    resid = y - X@beta
    # 常规 WLS 协方差（异方差稳健）；用 pinv + 截断负对角元避免数值 NaN
    meat = (X*w[:,None]*resid[:,None]).T @ (X*w[:,None]*resid[:,None])
    inv = np.linalg.pinv(X.T@(X*w[:,None]))
    cov = inv@meat@inv
    vd = np.clip(np.diag(cov), 0, None)
    return beta, np.sqrt(vd)

def rubin(ests, ses):
    Q=np.array(ests).T; U=np.array(ses).T**2; M=Q.shape[1]
    qbar=Q.mean(1); ubar=U.mean(1)
    b=((Q-qbar[:,None])**2).sum(1)/(M-1)
    T=ubar+(1+1/M)*b
    dfree=(M-1)*(1+ubar/((1+1/M)*b))**2
    dfree[~np.isfinite(dfree)|(b<=0)]=1e6
    t=qbar/np.sqrt(T); p=2*stats.t.sf(np.abs(t),dfree)
    return qbar,np.sqrt(T),t,p,dfree

# ---- 主分析 + BRR 重复权重方差 ----------------------------------------------
# Fay BRR: SE^2 = (1/(R*(1-rho)^2)) * sum (b_r - b_bar)^2, rho=0.5
rho=0.5; R=len(rw)
def brr_se(X, y, w_full, rw_mat):
    b_full,_ = wls(X,y,w_full)
    bs=[]
    for r in range(rw_mat.shape[1]):
        wr = ((1-rho)*w_full + rho*rw_mat[:,r])
        wr = wr/wr.mean()
        br,_ = wls(X,y,wr)
        bs.append(br)
    bs=np.array(bs)
    var = ((bs-bs.mean(0))**2).sum(0)/(R*(1-rho)**2)
    return np.sqrt(var)

results={}
for label in ["main","mod"]:
    X,names = design(label)
    ests=[]; ses=[]; brr_ses=[]
    for pvc in pv:
        y=d[pvc].astype(float).values
        b,se=wls(X,y,w_main)
        ests.append(b); ses.append(se)
        brr_ses.append(brr_se(X,y,w_main,rw_norm))
    qbar,se_r,t,p,dfree=rubin(ests,ses)
    # BRR 方差跨 PV 合并（用 Rubin 规则处理 BRR 得到的 SE）
    _,_,_,p_brr,_=rubin(ests,brr_ses)
    q2,se_brr,_,_,_=rubin(ests,brr_ses)
    res=pd.DataFrame({"term":names,"est":qbar,"se_model":se_r,
                      "se_brr":se_brr,"t":qbar/se_brr,"p_brr":p_brr})
    results[label]=res
    fn="model2_invertedU.csv" if label=="main" else "model3_moderation.csv"
    res.to_csv(f"{OUT}/{fn}",index=False,encoding="utf-8")
    r=res.copy()
    r["sig"]=np.where(r.p_brr<0.001,"***",np.where(r.p_brr<0.01,"**",np.where(r.p_brr<0.05,"*","")))
    print(f"\n===== {label}（5PV Rubin合并 + BRR SE，国家内归一化权重）=====")
    print(r.to_string(index=False,float_format=lambda x:f"{x:.4f}"),flush=True)

# 倒U判定 + 顶点位置诊断
r=results["main"]
b2=float(r.loc[r.term=="ai2","est"].iloc[0]); p2=float(r.loc[r.term=="ai2","p_brr"].iloc[0])
b1=float(r.loc[r.term=="ai","est"].iloc[0]); vertex_c=-b1/(2*b2)
# 顶点在原始 AIUSESCH 尺度 = 中心化顶点 + 均值；AIUSESCH 理论/实际范围约 [-1.6, 2.4]
am = df["AIUSESCH"].mean(); vertex_raw = vertex_c + am
lo, hi = df["AIUSESCH"].quantile(0.01), df["AIUSESCH"].quantile(0.99)
print("\n===== 倒U判定 =====")
print(f"ai2(二次项)系数={b2:.4f}  p(BRR)={p2:.4g}  {'显著' if p2<0.05 else '不显著'}")
print(f"顶点位置: 中心化尺度={vertex_c:.3f}  原始AIUSESCH尺度={vertex_raw:.3f}")
print(f"AIUSESCH 观测1-99%分位范围=[{lo:.3f},{hi:.3f}]  顶点是否落在数据范围内: {lo<=vertex_raw<=hi}")
if b2<0 and p2<0.05: print(f"[结论] 支持倒U型；顶点(原始尺度)={vertex_raw:.3f}")
elif b2>0 and p2<0.05:
    print("[结论] 显著【U型】(二次项为正)，与论文预期的倒U相反")
    if not (lo<=vertex_raw<=hi):
        print("[重要限定] 顶点落在观测数据范围之外，实际观测区间内呈单调下降，'U型'缺乏数据支撑，不可按倒U叙事解读")
else: print("[结论] 二次项不显著 => 无证据支持倒U，与论文预期不符")
# 保存倒U判定摘要
with open(f"{OUT}/verdict.txt","w",encoding="utf-8") as f:
    f.write(f"ai2系数={b2:.4f}; p(BRR)={p2:.4g}; 顶点原始尺度={vertex_raw:.3f}; "
            f"观测范围=[{lo:.3f},{hi:.3f}]; 顶点在范围内={lo<=vertex_raw<=hi}\n")

# 稳健性1：ST438四题均值
d2=d.copy()
d2["ai4"]=d2[["ST438Q01DA","ST438Q02DA","ST438Q03DA","ST438Q04DA"]].mean(axis=1,skipna=True)
d2=d2.dropna(subset=["ai4"]); d2["ai4c"]=d2.ai4-d2.ai4.mean(); d2["ai4c2"]=d2.ai4c**2
Xr=np.column_stack([np.ones(len(d2)),d2.ai4c,d2.ai4c2,d2.escs_c,d2.gender])
wr=norm_weights(d2.W_FSTUWT.values,d2.CNT.values)
ests=[wls(Xr,d2[pvc].astype(float).values,wr)[0] for pvc in pv]
q,se_r,_,p_r,_=rubin(ests,[wls(Xr,d2[pvc].astype(float).values,wr)[1] for pvc in pv])
print("\n===== 稳健性1：ST438四题均值替代AIUSESCH =====")
for nm,bb,pp in zip(["const","ai4","ai4^2","escs","gender"],q,p_r):
    print(f"  {nm:10s} {bb:8.4f}  p={pp:.4g}")
pd.DataFrame({"term":["const","ai4","ai4sq","escs","gender"],"est":q,"p":p_r}).to_csv(
    f"{OUT}/robust1_ai4mean.csv",index=False,encoding="utf-8")
print(f"\n总耗时(秒)={time.time()-t0:.1f}  结果目录={OUT}",flush=True)
