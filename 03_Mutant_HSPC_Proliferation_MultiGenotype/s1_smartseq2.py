import os, sys, warnings, numpy as np, pandas as pd, scanpy as sc
from scipy import stats
warnings.filterwarnings("ignore"); sc.settings.verbosity=0
sys.path.insert(0,os.path.dirname(__file__))
from genesets_mouse import S_PHASE,G2M_PHASE,PROLIFERATION_CORE
DATA=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); RES=os.path.join(DATA,"results"); os.makedirs(RES,exist_ok=True)
def proc(ad):
    sc.pp.filter_cells(ad,min_genes=500); sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
    s=[g for g in S_PHASE if g in ad.var_names]; g2=[g for g in G2M_PHASE if g in ad.var_names]
    sc.tl.score_genes_cell_cycle(ad,s_genes=s,g2m_genes=g2)
    pc=[g for g in PROLIFERATION_CORE if g in ad.var_names]; sc.tl.score_genes(ad,pc,score_name="prolif_core")
    return ad
plates={"WT":["GSE124822_LTHSC-WT_476_RNA_counts.txt.gz","GSE124822_LTHSC-WT_846_RNA_counts.txt.gz"],
        "Tet2":["GSE124822_LTHSC-TET2_30_RNA_counts.txt.gz","GSE124822_LTHSC_TET2_272_RNA_counts.txt.gz"],
        "Dnmt3a":["GSE124822_LTHSC-DNMT3A_675_RNA_counts.txt.gz"]}
rows=[]; obs={}
for geno,files in plates.items():
    mats=[pd.read_csv(os.path.join(DATA,f),sep="\t",index_col=0) for f in files]
    mats=[m.groupby(m.index).sum() for m in mats]
    genes=sorted(set().union(*[set(m.index) for m in mats]))
    expr=pd.concat([m.reindex(genes).fillna(0) for m in mats],axis=1)
    ad=sc.AnnData(expr.T.astype("float32")); ad.var_names_make_unique(); ad=proc(ad)
    obs[geno]=ad.obs.copy()
    ph=ad.obs["phase"].value_counts(normalize=True)*100
    rows.append(dict(dataset="GSE124822_SmartSeq2",group=geno,n_cells=ad.n_obs,
        pct_G1=ph.get("G1",0),pct_S=ph.get("S",0),pct_G2M=ph.get("G2M",0),
        pct_cycling=(ad.obs["phase"]!="G1").mean()*100,mean_prolif=ad.obs["prolif_core"].mean()))
pd.DataFrame(rows).to_csv(os.path.join(RES,"scrna_smartseq2_summary.csv"),index=False)
for g,o in obs.items(): o[["phase","S_score","G2M_score","prolif_core"]].to_csv(os.path.join(RES,f"obs_smartseq2_{g}.csv"))
print(pd.DataFrame(rows).to_string(index=False))
for mut in ["Tet2","Dnmt3a"]:
    w=obs["WT"]; m=obs[mut]
    tab=np.array([[(m["phase"]!="G1").sum(),(m["phase"]=="G1").sum()],[(w["phase"]!="G1").sum(),(w["phase"]=="G1").sum()]])
    chi2,p,_,_=stats.chi2_contingency(tab); u,pu=stats.mannwhitneyu(m["prolif_core"],w["prolif_core"])
    print(f"{mut} vs WT: cycling WT={ (w['phase']!='G1').mean()*100:.1f}% MUT={(m['phase']!='G1').mean()*100:.1f}% chi2_p={p:.3g} prolif_MWU_p={pu:.3g}")
