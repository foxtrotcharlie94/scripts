import os,numpy as np,pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
DATA=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); RES=os.path.join(DATA,"results")
obs={}
for k in ["smartseq2_WT","smartseq2_Tet2","smartseq2_Dnmt3a","tenx_WT","tenx_Jak2V617F"]:
    obs[k]=pd.read_csv(os.path.join(RES,f"obs_{k}.csv"),index_col=0)
def frac(o): 
    ph=o["phase"].value_counts(normalize=True)*100
    return ph.get("G1",0),ph.get("S",0),ph.get("G2M",0),(o["phase"]!="G1").mean()*100,o["prolif_core"].mean()
rows=[]; comp=[]
order=[("GSE124822 Smart-seq2 LT-HSC","smartseq2_WT","WT"),
       ("GSE124822 Smart-seq2 LT-HSC","smartseq2_Tet2","Tet2"),
       ("GSE124822 Smart-seq2 LT-HSC","smartseq2_Dnmt3a","Dnmt3a"),
       ("GSE263101 10x LSK","tenx_WT","WT"),
       ("GSE263101 10x LSK","tenx_Jak2V617F","Jak2V617F")]
for ds,k,g in order:
    g1,s,g2m,cyc,pr=frac(obs[k])
    rows.append(dict(dataset=ds,genotype=g,n_cells=len(obs[k]),pct_G1=round(g1,1),pct_S=round(s,1),
                     pct_G2M=round(g2m,1),pct_cycling=round(cyc,1),mean_prolif=round(pr,3)))
def cc(wk,mk,ds,mut):
    w,m=obs[wk],obs[mk]
    tab=np.array([[(m["phase"]!="G1").sum(),(m["phase"]=="G1").sum()],[(w["phase"]!="G1").sum(),(w["phase"]=="G1").sum()]])
    chi2,p,_,_=stats.chi2_contingency(tab); u,pu=stats.mannwhitneyu(m["prolif_core"],w["prolif_core"],alternative="two-sided")
    comp.append(dict(dataset=ds,comparison=f"{mut} vs WT",cycling_WT=round((w['phase']!='G1').mean()*100,1),
        cycling_MUT=round((m['phase']!='G1').mean()*100,1),delta_cycling=round(((m['phase']!='G1').mean()-(w['phase']!='G1').mean())*100,1),
        chi2_p=p,prolif_WT=round(w['prolif_core'].mean(),3),prolif_MUT=round(m['prolif_core'].mean(),3),MWU_p=pu))
cc("smartseq2_WT","smartseq2_Tet2","GSE124822","Tet2")
cc("smartseq2_WT","smartseq2_Dnmt3a","GSE124822","Dnmt3a")
cc("tenx_WT","tenx_Jak2V617F","GSE263101","Jak2V617F")
fr=pd.DataFrame(rows); cp=pd.DataFrame(comp)
fr.to_csv(os.path.join(RES,"MASTER_scrna_phase_fractions.csv"),index=False)
cp.to_csv(os.path.join(RES,"MASTER_scrna_comparisons.csv"),index=False)
print(fr.to_string(index=False)); print(); print(cp.to_string(index=False))

# ---- Figures ----
plt.rcParams.update({"font.size":10,"figure.dpi":130})
# Fig1: stacked phase fractions
fig,ax=plt.subplots(figsize=(8,4.2))
labels=[f"{r['genotype']}\n(n={r['n_cells']})" for _,r in fr.iterrows()]
G1=fr["pct_G1"];S=fr["pct_S"];G2M=fr["pct_G2M"]
ax.bar(labels,G1,color="#c9ced6",label="G1")
ax.bar(labels,S,bottom=G1,color="#4c9be8",label="S")
ax.bar(labels,G2M,bottom=G1+S,color="#e8744c",label="G2M")
for i,(_,r) in enumerate(fr.iterrows()): ax.text(i,101,f"{r['pct_cycling']:.0f}%",ha="center",fontsize=9,fontweight="bold")
ax.axvline(2.5,color="k",lw=0.8,ls="--"); ax.set_ylabel("% of cells"); ax.set_ylim(0,108)
ax.set_title("Cell-cycle phase by genotype (% above bars = cycling, S+G2M)")
ax.legend(ncol=3,loc="lower center",bbox_to_anchor=(0.5,-0.28),frameon=False)
ax.text(1,-16,"Smart-seq2 LT-HSC (GSE124822)",ha="center",fontsize=8,color="#555")
ax.text(3.5,-16,"10x LSK (GSE263101)",ha="center",fontsize=8,color="#555")
plt.tight_layout(); plt.savefig(os.path.join(RES,"fig1_phase_fractions.png"),bbox_inches="tight")
# Fig2: cycling fraction with deltas
fig,ax=plt.subplots(figsize=(7,4))
cols={"WT":"#7f8c8d","Tet2":"#2980b9","Dnmt3a":"#16a085","Jak2V617F":"#c0392b"}
ax.bar(range(len(fr)),fr["pct_cycling"],color=[cols[g] for g in fr["genotype"]])
ax.set_xticks(range(len(fr))); ax.set_xticklabels([r['genotype'] for _,r in fr.iterrows()],rotation=20)
ax.set_ylabel("% cycling (S+G2M)"); ax.set_title("Proliferating fraction by genotype")
for i,v in enumerate(fr["pct_cycling"]): ax.text(i,v+0.6,f"{v:.1f}",ha="center",fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(RES,"fig2_cycling_fraction.png"),bbox_inches="tight")
# Fig3: prolif score violins
fig,ax=plt.subplots(figsize=(7,4))
data=[obs[k]["prolif_core"].values for _,k,_ in order]; names=[g for _,_,g in order]
parts=ax.violinplot(data,showmeans=True,showextrema=False)
ax.set_xticks(range(1,len(names)+1)); ax.set_xticklabels(names,rotation=20)
ax.set_ylabel("Proliferation-core score (per cell)"); ax.set_title("Per-cell proliferation score by genotype")
ax.axvline(3.5,color="k",lw=0.8,ls="--")
plt.tight_layout(); plt.savefig(os.path.join(RES,"fig3_prolif_score.png"),bbox_inches="tight")
print("\nfigures + master tables written to results/")
