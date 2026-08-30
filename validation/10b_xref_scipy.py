"""Referência cruzada independente: scipy/numpy recalculam Benford e
uniformidade de dígitos sobre os MESMOS dados exportados pelo pacote."""
import csv, math, collections
import numpy as np
from scipy import stats
import os

here = os.path.dirname(os.path.abspath(__file__))
dados = collections.defaultdict(list)
with open(os.path.join(here, "xref", "dados.csv")) as f:
    for row in csv.DictReader(f):
        dados[row["dataset"]].append(int(row["valor"]))

def benford_exp(d):
    if d == 1:
        return np.array([math.log10(1 + 1/k) for k in range(1, 10)])
    return np.array([sum(math.log10(1 + 1/(10*k + j)) for k in range(1, 10))
                     for j in range(0, 10)])

def digits(xs, kind):
    out = []
    for x in xs:
        s = str(abs(x))
        if kind == "benford1" and x >= 1:   out.append(int(s[0]))
        elif kind == "benford2" and x >= 10: out.append(int(s[1]))
        elif kind == "last" and x >= 10:     out.append(int(s[-1]))
        elif kind == "penultimate" and x >= 100: out.append(int(s[-2]))
    return out

jl = {}
with open(os.path.join(here, "xref", "julia_resultados.csv")) as f:
    for row in csv.DictReader(f):
        jl[(row["dataset"], row["teste"])] = row

print("=" * 106)
print("REFERÊNCIA CRUZADA — ElectionForensics.jl (Julia) vs. scipy/numpy (Python)")
print("=" * 106)
hdr = f"{'dataset':<14}{'teste':<14}{'χ² julia':>14}{'χ² scipy':>14}{'p julia':>13}{'p scipy':>13}{'Δrel χ²':>11}{'Δ p':>11}"
print(hdr); print("-" * 106)

maxrel = 0.0
for ds in sorted(dados):
    xs = dados[ds]
    for teste, kind, labels in (("benford1","benford1",range(1,10)),
                                ("benford2","benford2",range(0,10)),
                                ("last","last",range(0,10)),
                                ("penultimate","penultimate",range(0,10))):
        d = digits(xs, kind)
        n = len(d)
        cnt = np.array([d.count(k) for k in labels], dtype=float)
        if teste == "benford1":   exp = benford_exp(1) * n
        elif teste == "benford2": exp = benford_exp(2) * n
        else:                     exp = np.full(len(cnt), n / 10.0)
        chi2 = float(np.sum((cnt - exp) ** 2 / exp))
        df = len(cnt) - 1
        p = float(stats.chi2.sf(chi2, df))
        # conferência dupla: scipy.stats.chisquare deve dar o mesmo χ²
        chi2b, pb = stats.chisquare(cnt, exp)
        assert abs(chi2 - chi2b) < 1e-8 * max(1, chi2), (chi2, chi2b)
        r = jl[(ds, teste)]
        cj, pj = float(r["chi2"]), float(r["pvalue"])
        rel = abs(chi2 - cj) / max(abs(cj), 1e-12)
        maxrel = max(maxrel, rel)
        print(f"{ds:<14}{teste:<14}{cj:>14.6f}{chi2:>14.6f}{pj:>13.4e}{p:>13.4e}{rel:>11.2e}{abs(p-pj):>11.2e}")

print("-" * 106)
print(f"discrepância relativa máxima em χ²: {maxrel:.3e}   (tolerância numérica ~1e-12)")

# MAD independente
print()
print("MAD (Nigrini) — conferência independente:")
for ds in sorted(dados):
    for teste, kind, labels, d in (("benford1","benford1",range(1,10),1),
                                   ("benford2","benford2",range(0,10),2)):
        dd = digits(dados[ds], kind); n = len(dd)
        obs = np.array([dd.count(k) for k in labels]) / n
        mad = float(np.mean(np.abs(obs - benford_exp(d))))
        mj = float(jl[(ds, teste)]["mad"])
        print(f"  {ds:<14}{teste:<12} julia={mj:.10f}  python={mad:.10f}  Δ={abs(mad-mj):.2e}")
