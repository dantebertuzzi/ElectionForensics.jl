#!/usr/bin/env python3
"""Baixa e normaliza os dados eleitorais reais usados na validação.

    python3 validation/fetch_dados_reais.py [dir_saida]

Duas fontes públicas de proveniência documentada:

(a) OpenElections — resultados por precinct das eleições gerais americanas de
    2020, transcritos dos boletins oficiais das secretarias eleitorais
    estaduais. https://github.com/openelections
    Pleitos sem alegação séria de fraude sistemática ⇒ testam CALIBRAÇÃO.

(b) dkobak/elections — protocolos das seções (UIK) das eleições federais
    russas, coletados por Sergey Shpilkin. É o dataset de Kobak, Shpilkin &
    Pshenichnikov (2016), "Integer percentages as electoral falsification
    fingerprints", Annals of Applied Statistics 10(1):54-73.
    Anomalia documentada e revisada por pares ⇒ testa PODER.

Saída: um CSV por pleito com (votos_do_vencedor, total_de_votos) por seção.

Os CSVs do OpenElections NÃO têm schema uniforme entre estados (AK não tem
`county`, RI traz `absentee_votes` numa coluna à parte, VT tem `town`), por
isso a leitura é por NOME de coluna e cada estado passa por uma auditoria
impressa: quais rótulos foram tratados como candidato, quais como linha
administrativa, e se a soma bate com a linha-resumo da própria fonte.
"""
import csv, collections, io, os, re, sys, urllib.request, zipfile

OUT = sys.argv[1] if len(sys.argv) > 1 else "validation/dados_reais"
# Estados cujo arquivo 2020 general/precinct existe no OpenElections.
UF = ["ak", "al", "az", "ca", "co", "ct", "fl", "hi", "ia", "id", "il", "in",
      "la", "ma", "me", "mi", "mn", "mo", "ms", "mt", "nd", "ne", "nh", "nm",
      "nv", "ny", "oh", "ok", "or", "pa", "ri", "sc", "sd", "tn", "tx", "ut",
      "vt", "wa", "wv", "wy"]
OE = ("https://raw.githubusercontent.com/openelections/openelections-data-{st}"
      "/master/2020/20201103__{st}__general__precinct.csv")
RU = "https://raw.githubusercontent.com/dkobak/elections/master/data/2012.csv.zip"

# Linhas administrativas: não são votos num candidato e entrariam em
# duplicidade no denominador. O padrão casa o rótulo INTEIRO — uma busca por
# substring removeria "Don Blankenship" (candidato real de 2020, presente em
# AK e VT) por conter "blank". Write-ins são votos de verdade e ficam.
ADMIN = re.compile(
    r"(total(s)?([ _-]*(votes|ballots)([ _-]*cast)?)?"
    r"|ballots([ _-]*cast)?|registered([ _-]*voters?)?"
    r"|blanks?([ _-]*votes?)?|over[ _-]*votes?|under[ _-]*votes?"
    r"|void|spoiled|turnout|exhausted([ _-]*ballots)?|scattering|scatter)", re.I)
# subconjunto que representa o total do precinct, usado só para conferência
RESUMO = re.compile(r"total(s)?([ _-]*(votes|ballots)([ _-]*cast)?)?", re.I)
# Brancos e nulos: não são votos em candidato, mas entram no total de cédulas.
# Servem para conferir a identidade `resumo = Σcandidatos + residual`.
RESIDUAL = re.compile(r"blanks?|spoiled|over[ _-]*votes?|under[ _-]*votes?", re.I)
# ME e outros estados trazem linhas de SUBTOTAL por condado disfarçadas de
# precinct ("ANDROSCOGGIN Total"). Incluí-las duplicaria seções inteiras.
AGREGADO = re.compile(r".*\b(total|subtotal|county[ _-]*wide|statewide)\b.*", re.I)


def baixar(url, cache=None):
    """Baixa com cache em disco: reexecutar a auditoria não repete a rede."""
    if cache and os.path.exists(cache):
        with open(cache, "rb") as f:
            return f.read()
    with urllib.request.urlopen(url, timeout=180) as r:
        dados = r.read()
    if cache:
        os.makedirs(os.path.dirname(cache), exist_ok=True)
        with open(cache, "wb") as f:
            f.write(dados)
    return dados


def inteiro(x):
    x = (x or "").strip().replace(",", "")
    try:
        return int(x)
    except ValueError:
        return None


def eua(st, destino):
    try:
        bruto = baixar(OE.format(st=st),
                       os.path.join(destino, "raw", f"{st}.csv")).decode("utf-8", "replace")
    except Exception as e:
        print(f"  {st.upper()}: falhou ({e})")
        return
    agg = collections.defaultdict(dict)
    tot, cands = collections.Counter(), collections.Counter()
    resumo, residual = collections.Counter(), collections.Counter()
    admins = set()
    for row in csv.DictReader(io.StringIO(bruto)):
        if (row.get("office") or "").strip() != "President":
            continue
        cand = (row.get("candidate") or "").strip()
        v = inteiro(row.get("votes"))
        if v is None:
            continue
        # RI reporta votos presenciais e por correspondência em colunas separadas
        v += inteiro(row.get("absentee_votes")) or 0
        chave = tuple((row.get(c) or "").strip() for c in ("county", "town", "precinct"))
        # descarta linhas agregadas que se apresentam como precinct
        if any(AGREGADO.fullmatch(x) for x in chave if x):
            continue
        if not cand or ADMIN.fullmatch(cand):
            admins.add(cand)
            if cand and RESUMO.fullmatch(cand):
                resumo[chave] += v
            elif cand and RESIDUAL.fullmatch(cand):
                residual[chave] += v
            continue
        agg[chave][cand] = agg[chave].get(cand, 0) + v
        tot[chave] += v
        cands[cand] += v
    if not cands:
        print(f"  {st.upper()}: sem dados de President")
        return
    venc = cands.most_common(1)[0][0]
    recs = [(agg[k].get(venc, 0), tot[k], k[0], k[2]) for k in agg if tot[k] > 0]
    soma = sum(cands.values())
    top2 = sum(v for _, v in cands.most_common(2)) / soma
    share_venc = cands[venc] / soma

    # CRIVO DE QUALIDADE. Vários arquivos do OpenElections trazem o mesmo
    # candidato com grafias diferentes entre condados ("BIDEN / HARRIS",
    # "Joseph R. Biden", "Biden, Joseph R."). Nesses casos os votos ficam
    # espalhados entre rótulos e `winner_votes` por precinct fica errado —
    # analisar isso produziria "anomalia" que é só erro de extração.
    # Numa presidencial americana os dois primeiros colocados somam ~97 %
    # e o vencedor estadual fica entre 40 % e 80 %.
    motivos = []
    if top2 < 0.85:
        motivos.append(f"top-2 = {100*top2:.0f} % (< 85 %: rótulos fragmentados)")
    if not (0.40 <= share_venc <= 0.80):
        motivos.append(f"vencedor com {100*share_venc:.0f} % (fora de 40–80 %)")
    if len(recs) < 150:
        motivos.append(f"apenas {len(recs)} precincts")

    print(f"  {st.upper()}: {len(recs):>5} precincts · {len(cands)} rótulos · "
          f"vencedor {venc} ({100 * share_venc:.1f} %) · top-2 {100*top2:.1f} %")
    if admins:
        print(f"        administrativos: {sorted(a for a in admins if a)}")
    # Integridade: onde a fonte traz a própria linha-resumo, deve valer
    # `resumo = Σcandidatos + brancos/nulos`. É a checagem de que a extração
    # não perdeu nem duplicou votos.
    comuns = [k for k in resumo if k in tot]
    if comuns:
        ok = sum(1 for k in comuns if resumo[k] == tot[k] + residual[k])
        if ok < len(comuns):
            motivos.append(f"integridade {ok}/{len(comuns)} contra a linha-resumo")
        print(f"        integridade (resumo = Σcandidatos + brancos/nulos): "
              f"{ok}/{len(comuns)}" + ("  ok" if ok == len(comuns) else "  ⚠"))

    if motivos:
        print(f"        ✗ EXCLUÍDO: {'; '.join(motivos)}")
        return (st, len(recs), len(cands), share_venc, top2, "excluido")

    with open(os.path.join(destino, f"us2020_{st}.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["winner_votes", "total_votes", "county", "precinct"])
        w.writerows(recs)
    return (st, len(recs), len(cands), share_venc, top2, "ok")


def russia(destino):
    # colunas do protocolo: 12 válidos, 22 Zyuganov, 25 Putin
    I_VALID, I_ZYUG, I_PUTIN = 12, 22, 25
    z = zipfile.ZipFile(io.BytesIO(baixar(RU)))
    nome = [n for n in z.namelist() if n.endswith(".csv")][0]
    r = csv.reader(io.TextIOWrapper(z.open(nome), encoding="utf-8", newline=""))
    hdr = next(r)
    assert "Путин" in hdr[I_PUTIN], f"layout inesperado: {hdr[I_PUTIN]!r}"
    recs = []
    for x in r:
        try:
            valid, putin, zyug = int(x[I_VALID]), int(x[I_PUTIN]), int(x[I_ZYUG])
        except (ValueError, IndexError):
            continue
        if valid > 0 and 0 <= putin <= valid and 0 <= zyug <= valid:
            recs.append((putin, zyug, valid, x[0]))
    with open(os.path.join(destino, "russia2012.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["putin", "zyuganov", "valid", "region"])
        w.writerows(recs)
    print(f"  Rússia 2012: {len(recs)} seções (UIK)")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    print("EUA 2020 (OpenElections):")
    qualidade = [eua(st, OUT) for st in UF]
    qualidade = [q for q in qualidade if q]
    with open(os.path.join(OUT, "qualidade.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["uf", "precincts", "rotulos", "share_vencedor", "top2", "status"])
        w.writerows(qualidade)
    apr = [q for q in qualidade if q[5] == "ok"]
    print(f"\n  aprovados no crivo de qualidade: {len(apr)}/{len(qualidade)} estados")
    print("Rússia 2012 (Shpilkin / Kobak et al.):")
    russia(OUT)
    print(f"\nPronto. Rode:  julia -t auto validation/16_dados_reais.jl {OUT}")
