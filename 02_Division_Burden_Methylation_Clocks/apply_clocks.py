import gzip, pickle, re
import pandas as pd
import numpy as np

d = pickle.load(open('/tmp/gse44117/clocks_mm9.pkl', 'rb'))
blood = d['blood']   # list of (chr, pos_mm9_1based-ish, weight)  -- pos is (lifted excel pos - 1)
wlmt = d['wlmt']

# Build lookup: for each clock cpg, the two candidate bed-start coordinates (offset -1 and 0)
def build_targets(cpglist):
    targets = {}  # (chr, cand_pos) -> list of (clock_index)
    for i, (c, pos, w) in enumerate(cpglist):
        for off in (-1, 0):
            targets.setdefault((c, pos + off), []).append(i)
    return targets

blood_targets = build_targets(blood)
wlmt_targets = build_targets(wlmt)

samples = {
    "HSC_young_1": "GSM1079935_RRBS_cpgMethylation_Mouse_blood_HSC_young_1.RRBS.bed.gz",
    "HSC_young_2": "GSM1079939_RRBS_cpgMethylation_Mouse_blood_HSC_young_2.RRBS.bed.gz",
    "HSC_old_3": "GSM1079926_RRBS_cpgMethylation_Mouse_blood_HSC_old_3.RRBS.bed.gz",
    "HSC_old_4": "GSM1079927_RRBS_cpgMethylation_Mouse_blood_HSC_old_4.RRBS.bed.gz",
    "HSC_young_10_reconst_1": "GSM1079936_RRBS_cpgMethylation_Mouse_blood_HSC_young_10_reconst_1.RRBS.bed.gz",
    "HSC_young_10_reconst_2": "GSM1079937_RRBS_cpgMethylation_Mouse_blood_HSC_young_10_reconst_2.RRBS.bed.gz",
    "HSC_young_10_reconst_3": "GSM1079938_RRBS_cpgMethylation_Mouse_blood_HSC_young_10_reconst_3.RRBS.bed.gz",
    "HSC_old_10_reconst_2": "GSM1079924_RRBS_cpgMethylation_Mouse_blood_HSC_old_10_reconst_2.RRBS.bed.gz",
    "HSC_old_10_reconst_3": "GSM1079925_RRBS_cpgMethylation_Mouse_blood_HSC_old_10_reconst_3.RRBS.bed.gz",
}

group_map = {
    "HSC_young_1": "young_baseline", "HSC_young_2": "young_baseline",
    "HSC_old_3": "old_baseline", "HSC_old_4": "old_baseline",
    "HSC_young_10_reconst_1": "young_10_reconst", "HSC_young_10_reconst_2": "young_10_reconst",
    "HSC_young_10_reconst_3": "young_10_reconst",
    "HSC_old_10_reconst_2": "old_10_reconst", "HSC_old_10_reconst_3": "old_10_reconst",
}

def parse_sample(path, targets):
    # returns dict clock_index -> (meth, total) using first match found (prefer offset -1 which had more hits)
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
                        found[idx] = (m, t)
    return found

rows = []
for sample, fname in samples.items():
    path = f"/tmp/gse44117/{fname}"
    blood_found = parse_sample(path, blood_targets)
    wlmt_found = parse_sample(path, wlmt_targets)

    for clockname, cpglist, found in [("Blood", blood, blood_found), ("WLMT", wlmt, wlmt_found)]:
        n_total = len(cpglist)
        n_cov = len(found)
        pos_vals, neg_vals = [], []
        weighted_sum = 0.0
        weight_cov_sum = 0.0
        for idx, (c, pos, w) in enumerate(cpglist):
            if idx in found:
                m, t = found[idx]
                pct = 100.0 * m / t
                weighted_sum += w * pct
                weight_cov_sum += abs(w)
                if w > 0:
                    pos_vals.append(pct)
                else:
                    neg_vals.append(pct)
        rows.append({
            "sample": sample,
            "group": group_map[sample],
            "clock": clockname,
            "n_cpgs_total": n_total,
            "n_cpgs_covered": n_cov,
            "pct_covered": 100 * n_cov / n_total,
            "mean_pct_meth_pos_weight_cpgs": np.mean(pos_vals) if pos_vals else np.nan,
            "mean_pct_meth_neg_weight_cpgs": np.mean(neg_vals) if neg_vals else np.nan,
            "weighted_sum_score": weighted_sum,
            "weighted_sum_score_normalized": weighted_sum / weight_cov_sum if weight_cov_sum else np.nan,
        })
    print("done", sample)

df = pd.DataFrame(rows)
df.to_csv('/tmp/gse44117/clock_cpg_analysis.csv', index=False)
print(df.to_string(index=False))
