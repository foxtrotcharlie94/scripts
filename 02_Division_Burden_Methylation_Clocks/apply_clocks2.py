import gzip, pickle
import pandas as pd
import numpy as np

d = pickle.load(open('/tmp/gse44117/clocks2_mm9.pkl', 'rb'))
yomt = d['yomt']
liver = d['liver']

def build_targets(cpglist):
    targets = {}
    for i, (c, pos, w) in enumerate(cpglist):
        for off in (-1, 0):
            targets.setdefault((c, pos + off), []).append(i)
    return targets

yomt_targets = build_targets(yomt)
liver_targets = build_targets(liver)

samples = {
    "young_1": "GSM1079935_RRBS_cpgMethylation_Mouse_blood_HSC_young_1.RRBS.bed.gz",
    "young_2": "GSM1079939_RRBS_cpgMethylation_Mouse_blood_HSC_young_2.RRBS.bed.gz",
    "old_3": "GSM1079926_RRBS_cpgMethylation_Mouse_blood_HSC_old_3.RRBS.bed.gz",
    "old_4": "GSM1079927_RRBS_cpgMethylation_Mouse_blood_HSC_old_4.RRBS.bed.gz",
    "reconst_1": "GSM1079936_RRBS_cpgMethylation_Mouse_blood_HSC_young_10_reconst_1.RRBS.bed.gz",
    "reconst_2": "GSM1079937_RRBS_cpgMethylation_Mouse_blood_HSC_young_10_reconst_2.RRBS.bed.gz",
    "reconst_3": "GSM1079938_RRBS_cpgMethylation_Mouse_blood_HSC_young_10_reconst_3.RRBS.bed.gz",
}
young_samples = ["young_1", "young_2"]
old_samples = ["old_3", "old_4"]
reconst_samples = ["reconst_1", "reconst_2", "reconst_3"]

def parse_sample(path, targets):
    found = {}
    with gzip.open(path, 'rt') as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            key = (p[0], int(p[1]))
            if key in targets:
                m_str = p[3].strip("'")
                m, t = m_str.split('/')
                m, t = int(m), int(t)
                for idx in targets[key]:
                    if idx not in found:
                        found[idx] = 100.0 * m / t
    return found

results = {"YOMT": {}, "Liver": {}}
for sample, fname in samples.items():
    path = f"/tmp/gse44117/{fname}"
    results["YOMT"][sample] = parse_sample(path, yomt_targets)
    results["Liver"][sample] = parse_sample(path, liver_targets)
    print("parsed", sample)

all_rows = []
for clockname, cpglist in [("YOMT", yomt), ("Liver", liver)]:
    for idx, (c, pos, w) in enumerate(cpglist):
        vals = {}
        ok = True
        for s in samples:
            v = results[clockname][s].get(idx)
            if v is None:
                ok = False
                break
            vals[s] = v
        if not ok:
            continue
        young_mean = np.mean([vals[s] for s in young_samples])
        old_mean = np.mean([vals[s] for s in old_samples])
        reconst_mean = np.mean([vals[s] for s in reconst_samples])
        all_rows.append({
            "clock": clockname, "chr": c, "pos": pos, "weight": w,
            "delta_age_old_minus_young": old_mean - young_mean,
            "delta_reconst_minus_young": reconst_mean - young_mean,
        })

df2 = pd.DataFrame(all_rows)
df2.to_csv('/tmp/gse44117/cpg_deltas_clocks2.csv', index=False)
print(f"YOMT n={sum(df2.clock=='YOMT')}, Liver n={sum(df2.clock=='Liver')}")
