import os,sys,gzip,tarfile,warnings,numpy as np,pandas as pd
warnings.filterwarnings("ignore"); from scipy import stats
sys.path.insert(0,os.getcwd()); from genesets_mouse import S_PHASE,G2M_PHASE,PROLIFERATION_CORE
DATA=os.path.dirname(os.getcwd()); RES=sys.argv[1]; which=sys.argv[2]
PROL=sorted(set(g.lower() for g in S_PHASE+G2M_PHASE+PROLIFERATION_CORE))
KEY=["Mki67","Top2a","Pcna","Ccnb1","Ccnb2","Ccna2","Cdk1","Mcm2","Mcm5","Foxm1","Bub1","Aurka","Plk1","Rrm2"]
def deseq(counts,wt,mut):
    from pydeseq2.dds import DeseqDataSet; from pydeseq2.ds import DeseqStats
    counts=counts.groupby(counts.index).sum(); counts=counts[counts[wt+mut].sum(1)>=10]
    X=counts[wt+mut].T.round().astype(int)
    meta=pd.DataFrame({"condition":["WT"]*len(wt)+["MUT"]*len(mut)},index=wt+mut)
    dds=DeseqDataSet(counts=X,metadata=meta,design_factors="condition",quiet=True); dds.deseq2()
    st=DeseqStats(dds,contrast=["condition","MUT","WT"],quiet=True); st.summary()
    r=st.results_df; r.index=list(counts.index); return r
def report(r,name):
    r=r.dropna(subset=["log2FoldChange"]); low={g.lower():g for g in r.index}
    pg=[low[g] for g in PROL if g in low]; lp=r.loc[pg,"log2FoldChange"]
    U,pmw=stats.mannwhitneyu(lp,r["log2FoldChange"],alternative="greater"); W,pw=stats.wilcoxon(lp)
    up=(lp>0).sum()
    print(f"##### {name} #####  genes={len(r)} prolif={len(pg)}")
    print(f"  median LFC prolif={lp.median():+.3f} genome={r['log2FoldChange'].median():+.3f}")
    print(f"  competitive(up vs bg) p={pmw:.3g}  signedWilcoxon p={pw:.3g}  prolif up={up}/{len(pg)} ({100*up/len(pg):.0f}%)")
    print("  key:", ", ".join(f"{k} {r.loc[k,'log2FoldChange']:+.2f}/p{r.loc[k,'padj']:.1g}" for k in KEY if k in r.index))
    pd.DataFrame([dict(dataset=name,med_LFC_prolif=round(float(lp.median()),3),med_LFC_all=round(float(r['log2FoldChange'].median()),3),
        competitive_p=float(pmw),wilcoxon_p=float(pw),frac_up=round(float(up/len(pg)),2))]).to_csv(os.path.join(RES,f"rigor_{name.split()[0]}.csv"),index=False)

if which=="asxl1":
    ax=pd.read_csv(os.path.join(DATA,"GSE158184_Bulk_ASXL_counts.filt.txt.gz"),sep="\t",index_col=0)
    report(deseq(ax,["W2","W3","W4","W5"],["M2","M3","M4","M5"]),"Asxl1")
elif which=="jak2":
    feat={}
    with tarfile.open(os.path.join(DATA,"GSE263101_RAW.tar")) as tf:
        for line in gzip.open(tf.extractfile("GSM8185789_Vav1_iCre_features.tsv.gz"),"rt"):
            e,s=line.rstrip("\n").split("\t")[:2]; feat[e]=s
    jmap={"GSM3502680":"WT1","GSM3502681":"WT2","GSM3502682":"WT3","GSM3502683":"VF1","GSM3502684":"VF2","GSM3502685":"VF3"}
    cols={}
    with tarfile.open(os.path.join(DATA,"GSE123401_RAW.tar")) as tf:
        for m in tf.getmembers():
            g=m.name.split("_")[0]
            if g in jmap:
                d=pd.read_csv(gzip.open(tf.extractfile(m),"rt"),sep="\t",comment="#")
                c=d.iloc[:,[0,-1]]; c.columns=["gid","ct"]; c["sym"]=c["gid"].str.split(".").str[0].map(feat)
                cols[jmap[g]]=c.dropna(subset=["sym"]).groupby("sym")["ct"].sum()
    report(deseq(pd.DataFrame(cols),["WT1","WT2","WT3"],["VF1","VF2","VF3"]),"Jak2V617F")
