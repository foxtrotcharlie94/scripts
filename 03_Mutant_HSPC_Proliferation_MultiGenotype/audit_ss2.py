import os,sys,warnings,numpy as np,pandas as pd,scanpy as sc
warnings.filterwarnings("ignore"); sc.settings.verbosity=0
import statsmodels.formula.api as smf
sys.path.insert(0,os.getcwd()); from genesets_mouse import S_PHASE,G2M_PHASE,PROLIFERATION_CORE
DATA=os.path.dirname(os.getcwd()); RES=sys.argv[1]
plates={"WT_476":"GSE124822_LTHSC-WT_476_RNA_counts.txt.gz","WT_846":"GSE124822_LTHSC-WT_846_RNA_counts.txt.gz",
        "Tet2_30":"GSE124822_LTHSC-TET2_30_RNA_counts.txt.gz","Tet2_272":"GSE124822_LTHSC_TET2_272_RNA_counts.txt.gz",
        "Dnmt3a_675":"GSE124822_LTHSC-DNMT3A_675_RNA_counts.txt.gz"}
parts=[]
for plate,f in plates.items():
    d=pd.read_csv(os.path.join(DATA,f),sep="\t",index_col=0); d=d.groupby(d.index).sum()
    a=sc.AnnData(d.T.astype("float32")); a.var_names_make_unique()
    a.obs["n_genes"]=(a.X>0).sum(1); a.obs["total"]=a.X.sum(1)
    sc.pp.filter_cells(a,min_genes=500)
    sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
    s=[g for g in S_PHASE if g in a.var_names]; g2=[g for g in G2M_PHASE if g in a.var_names]
    sc.tl.score_genes_cell_cycle(a,s_genes=s,g2m_genes=g2)
    sc.tl.score_genes(a,[g for g in PROLIFERATION_CORE if g in a.var_names],score_name="prolif")
    o=a.obs.copy(); o["plate"]=plate; o["genotype"]=plate.split("_")[0]; parts.append(o)
df=pd.concat(parts,ignore_index=True); df["cycling"]=(df["phase"]!="G1").astype(int)
# per-plate table
tab=df.groupby(["genotype","plate"]).agg(n=("cycling","size"),med_genes=("n_genes","median"),
     med_total=("total","median"),cycling_pct=("cycling",lambda x:round(100*x.mean(),1)),
     S=("S_score","mean"),G2M=("G2M_score","mean"),prolif=("prolif","mean")).round(3)
print("=== per-plate (batch) QC + proliferation ==="); print(tab.to_string())
# depth confound: does detection drive cycling within WT?
wt=df[df.genotype=="WT"]
r=np.corrcoef(wt["n_genes"],wt["S_score"]+wt["G2M_score"])[0,1]
print(f"\nWithin-WT corr(n_genes, S+G2M score) = {r:.2f}")
print("median genes  WT=%.0f Tet2=%.0f Dnmt3a=%.0f"%(
  df[df.genotype=='WT'].n_genes.median(),df[df.genotype=='Tet2'].n_genes.median(),df[df.genotype=='Dnmt3a'].n_genes.median()))
# depth- and batch-controlled test: cycling ~ genotype + n_genes (Tet2/Dnmt3a vs WT)
df["ng_z"]=(df["n_genes"]-df["n_genes"].mean())/df["n_genes"].std()
for mut in ["Tet2","Dnmt3a"]:
    sub=df[df.genotype.isin(["WT",mut])].copy(); sub["mut"]=(sub.genotype==mut).astype(int)
    m1=smf.logit("cycling ~ mut",data=sub).fit(disp=0)
    m2=smf.logit("cycling ~ mut + ng_z",data=sub).fit(disp=0)
    print(f"\n{mut} vs WT:  unadjusted OR={np.exp(m1.params['mut']):.2f} p={m1.pvalues['mut']:.3g} | depth-adjusted OR={np.exp(m2.params['mut']):.2f} p={m2.pvalues['mut']:.3g}")
df.to_csv(os.path.join(RES,"ss2_audit_percell.csv"),index=False)
tab.to_csv(os.path.join(RES,"ss2_audit_perplate.csv"))
