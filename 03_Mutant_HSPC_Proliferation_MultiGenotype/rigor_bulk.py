import os,sys,gzip,tarfile,warnings,numpy as np,pandas as pd
warnings.filterwarnings("ignore")
from scipy import stats
sys.path.insert(0,os.getcwd()); from genesets_mouse import S_PHASE,G2M_PHASE,PROLIFERATION_CORE
DATA=os.path.dirname(os.getcwd()); RES=sys.argv[1]
PROL=sorted(set([g.lower() for g in S_PHASE+G2M_PHASE+PROLIFERATION_CORE]))
KEY=["Mki67","Top2a","Pcna","Ccnb1","Ccnb2","Ccna2","Cdk1","Mcm2","Mcm5","Foxm1","Bub1","Aurka","Plk1","Rrm2"]

def deseq(counts,wt,mut,name):
    from pydeseq2.dds import DeseqDataSet; from pydeseq2.ds import DeseqStats
    counts=counts.groupby(counts.index).sum()
    counts=counts[(counts[wt+mut].sum(1)>=10)]
    X=counts[wt+mut].T.round().astype(int)
    meta=pd.DataFrame({"condition":["WT"]*len(wt)+["MUT"]*len(mut)},index=wt+mut)
    dds=DeseqDataSet(counts=X,metadata=meta,design_factors="condition",quiet=True)
    dds.deseq2(); st=DeseqStats(dds,contrast=["condition","MUT","WT"],quiet=True); st.summary()
    r=st.results_df; r.index=[i for i in counts.index]
    return r

def report(r,name):
    r=r.dropna(subset=["log2FoldChange"])
    low={g.lower():g for g in r.index}
    pg=[low[g] for g in PROL if g in low]
    lfc_p=r.loc[pg,"log2FoldChange"]; lfc_all=r["log2FoldChange"]
    U,pmw=stats.mannwhitneyu(lfc_p,lfc_all,alternative="greater")  # prolif up vs background?
    W,pw=stats.wilcoxon(lfc_p)  # prolif LFCs != 0?
    print(f"\n##### {name} #####  (genes tested={len(r)}, prolif genes={len(pg)})")
    print(f"  median LFC: proliferation={lfc_p.median():+.3f}  genome={lfc_all.median():+.3f}")
    print(f"  competitive test (prolif up vs background) p={pmw:.3g}  | signed Wilcoxon p={pw:.3g}")
    nz=r.loc[pg]; up=(nz['log2FoldChange']>0).sum()
    print(f"  proliferation genes up:{up}/{len(pg)}  ({100*up/len(pg):.0f}%)")
    kk=[k for k in KEY if k in r.index]
    print("  key genes (LFC, padj):", ", ".join(f"{k} {r.loc[k,'log2FoldChange']:+.2f}/{r.loc[k,'padj']:.2g}" for k in kk))
    return dict(dataset=name,med_LFC_prolif=round(float(lfc_p.median()),3),med_LFC_all=round(float(lfc_all.median()),3),
                competitive_p=float(pmw),wilcoxon_p=float(pw),frac_up=round(float(up/len(pg)),2))

out=[]
# ASXL1 counts 4v4 (best replication)
ax=pd.read_csv(os.path.join(DATA,"GSE158184_Bulk_ASXL_counts.filt.txt.gz"),sep="\t",index_col=0)
out.append(report(deseq(ax,["W2","W3","W4","W5"],["M2","M3","M4","M5"],"Asxl1"),"Asxl1 (CD133+ 4v4, DESeq2)"))
# JAK2 counts 3v3
feat={}
with tarfile.open(os.path.join(DATA,"GSE263101_RAW.tar")) as tf:
    for line in gzip.open(tf.extractfile("GSM8185789_Vav1_iCre_features.tsv.gz"),"rt"):
        e,s=line.rstrip("\n").split("\t")[:2]; feat[e]=s
jmap={"GSM3502680":"WT1","GSM3502681":"WT2","GSM3502682":"WT3","GSM3502683":"VF1","GSM3502684":"VF2","GSM3502685":"VF3"}
cols={}
with tarfile.open(os.path.join(DATA,"GSE123401_RAW.tar")) as tf:
    for m in tf.getmembers():
        gsm=m.name.split("_")[0]
        if gsm in jmap:
            d=pd.read_csv(gzip.open(tf.extractfile(m),"rt"),sep="\t",comment="#")
            c=d.iloc[:,[0,-1]]; c.columns=["gid","ct"]; c["sym"]=c["gid"].str.split(".").str[0].map(feat)
            cols[jmap[gsm]]=c.dropna(subset=["sym"]).groupby("sym")["ct"].sum()
jak=pd.DataFrame(cols)
out.append(report(deseq(jak,["WT1","WT2","WT3"],["VF1","VF2","VF3"],"Jak2"),"Jak2V617F (LSK 3v3, DESeq2)"))
# TET2 normalized 3v3 -> per-gene Welch t + competitive on t-stats
t2=pd.read_csv(os.path.join(DATA,"TET2GSE96758_normalizedData.txt.gz"),sep="\t",index_col=0)
t2=t2.groupby(t2.index).mean()
wt=["WT1","WT2","WT3"]; mut=["T1","T2","T3"]
t,p=stats.ttest_ind(t2[mut],t2[wt],axis=1,equal_var=False)
res=pd.DataFrame({"t":t,"diff":t2[mut].mean(1)-t2[wt].mean(1)},index=t2.index).dropna()
low={g.lower():g for g in res.index}; pg=[low[g] for g in PROL if g in low]
U,pmw=stats.mannwhitneyu(res.loc[pg,"t"],res["t"],alternative="greater")
print(f"\n##### Tet2 (GSE96758 normalized 3v3) #####")
print(f"  median t: proliferation={res.loc[pg,'t'].median():+.2f} genome={res['t'].median():+.2f}  competitive p={pmw:.3g}")
kk=[k for k in KEY if k in res.index]
print("  key genes (mut-WT diff, t):", ", ".join(f"{k} {res.loc[k,'diff']:+.2f}/{res.loc[k,'t']:+.1f}" for k in kk))
out.append(dict(dataset="Tet2 (norm 3v3, t-test)",med_LFC_prolif=round(float(res.loc[pg,'t'].median()),3),
                med_LFC_all=round(float(res['t'].median()),3),competitive_p=float(pmw),wilcoxon_p=np.nan,frac_up=round(float((res.loc[pg,'diff']>0).mean()),2)))
pd.DataFrame(out).to_csv(os.path.join(RES,"bulk_proliferation_RIGOROUS.csv"),index=False)
print("\nsaved bulk_proliferation_RIGOROUS.csv")
