#!/usr/bin/env python3
"""Os concorrentes têm a mesma ressalva?

Roda as bibliotecas existentes de Benford, COM SEUS DEFAULTS, sobre as mesmas
29 eleições americanas legítimas de 2020 usadas em `16_dados_reais.jl`, e conta
quantas cada uma acusa. Sob calibração, o esperado é ~5 % a α = 0,05.

Bibliotecas executadas de fato:
  · benford_py     (PyPI, instalada)
  · benfordslaw    (PyPI, instalada)

Estatísticas de BenfordTests (CRAN) reimplementadas a partir do código-fonte R
(`Benford_tests.R`), porque não há R neste ambiente. A reimplementação é
conferida contra a definição no fonte, não contra a execução do pacote.

    python3 validation/20_concorrentes.py <dir_dos_csvs>
"""
import csv, glob, math, os, sys, warnings
import numpy as np
from scipy import stats

warnings.filterwarnings("ignore")
D = sys.argv[1] if len(sys.argv) > 1 else "validation/dados_reais"
ALPHA = 0.05


def pbenf(digits):
    """Distribuição de Benford dos primeiros `digits` dígitos (BenfordTests::pbenf)."""
    lo, hi = 10 ** (digits - 1), 10 ** digits
    return np.array([math.log10(1 + 1 / d) for d in range(lo, hi)])


def signifd(x, digits=1):
    """Primeiros `digits` dígitos significativos (BenfordTests::signifd)."""
    x = np.abs(np.asarray(x, dtype=float))
    x = x[x > 0]
    return np.floor(x / 10 ** np.floor(np.log10(x)) * 10 ** (digits - 1)).astype(int)


def counts(fd, digits):
    lo, hi = 10 ** (digits - 1), 10 ** digits
    return np.array([np.sum(fd == d) for d in range(lo, hi)], dtype=float)


# ── BenfordTests (CRAN), reimplementado do fonte ─────────────────────
def chisq_benftest(x, digits=1):
    fd = signifd(x, digits); n = len(fd)
    obs = counts(fd, digits) / n
    p0 = pbenf(digits)
    chi = n * np.sum((obs - p0) ** 2 / p0)
    return float(stats.chi2.sf(chi, len(p0) - 1))


def ks_benftest(x, digits=1):
    """Kolmogorov-Smirnov contra Benford. p-valor pela distribuição
    assintótica de Kolmogorov (o pacote R simula; usamos a assintótica)."""
    fd = signifd(x, digits); n = len(fd)
    obs = np.cumsum(counts(fd, digits) / n)
    exp_ = np.cumsum(pbenf(digits))
    d = np.max(np.abs(obs - exp_))
    lam = (math.sqrt(n) + 0.12 + 0.11 / math.sqrt(n)) * d
    p = 2 * sum((-1) ** (k - 1) * math.exp(-2 * k * k * lam * lam) for k in range(1, 100))
    return float(min(max(p, 0.0), 1.0))


def meandigit_benftest(x, digits=1):
    """Judge-Schechter, versão assintótica (BenfordTests::meandigit.benftest)."""
    fd = signifd(x, digits); n = len(fd)
    lo, hi = 10 ** (digits - 1), 10 ** digits
    seq = np.arange(lo, hi, dtype=float)
    p0 = pbenf(digits)
    mu_bed = float(seq @ p0)
    var_bed = float(((seq - mu_bed) ** 2) @ p0)
    mu_emp = float(np.mean(fd))
    a_star = abs(mu_emp - mu_bed) / (seq.max() - mu_bed)
    sd = math.sqrt(var_bed / n) / (seq.max() - mu_bed)
    return float(2 * (1 - stats.norm.cdf(a_star, 0, sd)))


# ── benford.analysis (CRAN): MAD + conformidade de Nigrini ───────────
# O pacote não devolve p-valor; o veredito é a faixa de MAD. Default: 2 dígitos.
MAD_NIGRINI_2D = (0.0012, 0.0018, 0.0022)   # 1º-2º dígitos (Nigrini 2012)


def benford_analysis_verdict(x, digits=2):
    fd = signifd(x, digits); n = len(fd)
    obs = counts(fd, digits) / n
    mad = float(np.mean(np.abs(obs - pbenf(digits))))
    t1, t2, t3 = MAD_NIGRINI_2D
    faixa = ("close" if mad <= t1 else "acceptable" if mad <= t2
             else "marginal" if mad <= t3 else "NONCONFORMITY")
    return mad, faixa


def ler(path):
    v = []
    with open(path) as f:
        for i, row in enumerate(csv.reader(f)):
            if i == 0 or len(row) < 2:
                continue
            try:
                a, b = int(row[0]), int(row[1])
            except ValueError:
                continue
            if b > 0 and 0 <= a <= b:
                v.append(a)
    return np.array(v)


def main():
    import benford as bf
    from benfordslaw import benfordslaw

    arquivos = sorted(glob.glob(os.path.join(D, "us2020_*.csv")))
    if not arquivos:
        sys.exit(f"nenhum us2020_*.csv em {D} — rode fetch_dados_reais.py antes")

    linhas = []
    for f in arquivos:
        uf = os.path.basename(f)[7:9].upper()
        v = ler(f)
        vv = v[v >= 10]
        if len(vv) < 150:
            continue
        r = {"uf": uf, "n": len(vv)}
        r["BenfordTests::chisq (1BL)"] = chisq_benftest(vv, 1)
        r["BenfordTests::ks (1BL)"] = ks_benftest(vv, 1)
        r["BenfordTests::meandigit (1BL)"] = meandigit_benftest(vv, 1)
        r["BenfordTests::chisq (2 díg.)"] = chisq_benftest(vv, 2)
        # benfordslaw, com seu default (1º dígito, teste qui-quadrado)
        bl = benfordslaw(alpha=ALPHA, method="chi2", verbose=0)
        r["benfordslaw (default)"] = float(bl.fit(vv)["P"])
        # benford.analysis: veredito por MAD, não p-valor
        mad, faixa = benford_analysis_verdict(vv, 2)
        r["benford.analysis (MAD 2 díg.)"] = faixa
        # benford_py: MAD do 1º dígito + faixa de conformidade de Nigrini
        t = bf.first_digits(vv.astype(float), digs=1, decimals=0,
                            confidence=None, verbose=False)
        mad1 = float(np.mean(np.abs(t["Found"].values - t["Expected"].values)))
        r["benford_py (MAD 1BL)"] = ("close" if mad1 <= 0.006 else
                                     "acceptable" if mad1 <= 0.012 else
                                     "marginal" if mad1 <= 0.015 else "NONCONFORMITY")
        linhas.append(r)

    testes_p = [k for k in linhas[0] if k not in ("uf", "n") and
                not isinstance(linhas[0][k], str)]
    testes_v = [k for k in linhas[0] if k != "uf" and isinstance(linhas[0][k], str)]

    print("=" * 108)
    print(f"OS CONCORRENTES ACUSAM ELEIÇÕES LEGÍTIMAS?  {len(linhas)} pleitos "
          f"americanos de 2020, defaults de cada biblioteca")
    print("=" * 108)
    largura = max(len(t) for t in testes_p + testes_v) + 2
    print(f"{'teste':<{largura}}{'rejeições @0,05':>18}{'taxa':>10}   {'KS vs U(0,1)':>14}")
    print("-" * 108)
    for t in testes_p:
        ps = np.array([l[t] for l in linhas])
        ks = float(stats.kstest(ps, "uniform").pvalue)
        rej = int(np.sum(ps < ALPHA))
        print(f"{t:<{largura}}{f'{rej}/{len(ps)}':>18}{rej/len(ps):>10.3f}   {ks:>14.2e}")
    print("-" * 108)
    for t in testes_v:
        naoconf = sum(1 for l in linhas if l[t] in ("NONCONFORMITY", "marginal"))
        print(f"{t:<{largura}}{f'{naoconf}/{len(linhas)}':>18}"
              f"{naoconf/len(linhas):>10.3f}   {'(sem p-valor)':>14}")
    print()
    print("Para comparação, medido em 16_dados_reais.jl sobre os mesmos dados:")
    print(f"{'ElectionForensics.jl 2BL null clássico':<{largura}}{'8/29':>18}{0.276:>10.3f}"
          f"{1.0e-4:>17.2e}")
    print(f"{'ElectionForensics.jl 2BL null reamostrado':<{largura}}{'0/29':>18}{0.0:>10.3f}"
          f"{0.16:>17.2e}")


if __name__ == "__main__":
    main()
