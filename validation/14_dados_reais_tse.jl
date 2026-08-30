# Validação em dados reais do TSE.
# NÃO EXECUTADO nesta revisão: o CDN do TSE respondeu HTTP 403 a clientes
# não-browser a partir do ambiente de auditoria. Rode localmente.
#
#   julia -t auto validation/14_dados_reais_tse.jl <caminho/votacao_secao_2022_XX.csv>
#
# Baixe em https://dadosabertos.tse.jus.br/dataset/resultados-2022
# (arquivo "votação por seção eleitoral"), descompacte e passe o CSV.
using ElectionForensics, Printf, Random

const CSV_PATH = length(ARGS) ≥ 1 ? ARGS[1] :
    error("uso: julia validation/14_dados_reais_tse.jl <votacao_secao_*.csv>")

# leitor mínimo do layout do TSE (latin-1, separador ';', campos entre aspas)
function ler_tse(path; cargo = "Presidente", turno = 2)
    linhas = eachline(open(path; lock = false))
    hdr = split(strip(first(Iterators.take(linhas, 1))), ';')
    error("Implemente a leitura conforme o layout vigente; veja o dicionário " *
          "de dados do TSE. Campos necessários: NR_TURNO, DS_CARGO, " *
          "SG_UF, NR_ZONA, NR_SECAO, NR_VOTAVEL, QT_VOTOS.")
end

# Esqueleto do que rodar depois de montar `votes` (votos do candidato por seção)
# e `totals` (votos válidos por seção):
#
#   res = forensics_report(votes, totals; B = 999, rng = MersenneTwister(1))
#
#   # Checagem de reprodutibilidade: dois RNGs distintos devem dar p-valores
#   # dentro do erro Monte Carlo (~0.014 com B = 999).
#   ps = [rozenas_test(votes, totals; B = 999, rng = MersenneTwister(k)).pvalue
#         for k in 1:10]
#   @printf("rozenas p: média %.4f, sd %.4f (esperado sd ≈ 0.014)\n",
#           mean(ps), std(ps))
#
#   # Diagnóstico de aplicabilidade dos nulls de dígito:
#   oom = log10(quantile(Float64.(votes), 0.95) / quantile(Float64.(votes), 0.05))
#   @printf("amplitude das contagens: %.2f décadas — 2BL clássica só é\n", oom)
#   @printf("defensável com tamanhos de seção log-normais dispersos.\n")
