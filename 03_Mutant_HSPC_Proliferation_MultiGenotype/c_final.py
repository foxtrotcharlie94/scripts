import os,sys,numpy as np,pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
RES=sys.argv[1]
# Final matched-compartment comparison (from computed results)
rows=[
 dict(genotype="Tet2",compartment="LT-HSC (Smart-seq2)",pop="LT-HSC",WT=42.4,MUT=28.8,p=0.0076,dir="down"),
 dict(genotype="Dnmt3a",compartment="LT-HSC (Smart-seq2)",pop="LT-HSC",WT=42.4,MUT=45.3,p=0.74,dir="ns"),
 dict(genotype="Jak2V617F",compartment="LT-HSC+MPP1 (10x)",pop="LT-HSC+MPP1",WT=45.7,MUT=52.9,p=5.5e-19,dir="up"),
]
df=pd.DataFrame(rows); df.to_csv(os.path.join(RES,"MASTER_matched_compartment.csv"),index=False)
fig,ax=plt.subplots(figsize=(8,4.6))
x=np.arange(len(df)); w=0.36
ax.bar(x-w/2,df["WT"],w,label="WT/control",color="#95a5a6")
cols={"down":"#2980b9","ns":"#bdc3c7","up":"#c0392b"}
ax.bar(x+w/2,df["MUT"],w,label="mutant",color=[cols[d] for d in df["dir"]])
for i,r in df.iterrows():
    ax.text(i-w/2,r.WT+0.6,f"{r.WT:.1f}",ha="center",fontsize=8)
    ax.text(i+w/2,r.MUT+0.6,f"{r.MUT:.1f}",ha="center",fontsize=8)
    star = "***" if r.p<1e-3 else ("**" if r.p<0.01 else ("*" if r.p<0.05 else "ns"))
    ax.text(i,max(r.WT,r.MUT)+3,star,ha="center",fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels([f"{r.genotype}\n{r.compartment}" for _,r in df.iterrows()],fontsize=9)
ax.set_ylabel("% cycling (S+G2M)"); ax.set_ylim(0,62)
ax.set_title("Proliferating fraction, mutant vs control — matched primitive compartment")
ax.legend(); plt.tight_layout(); plt.savefig(os.path.join(RES,"fig6_summary.png"),bbox_inches="tight",dpi=130)
print(df.to_string(index=False)); print("fig6 saved")
