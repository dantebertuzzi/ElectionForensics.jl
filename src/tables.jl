# ── integração com Tables.jl ─────────────────────────────────────────
#
# Cada resultado expõe a tabela que o usuário de fato quer plotar ou exportar:
# uma linha por dígito (Benford, último dígito) ou por fração (Rozenas). Os
# escalares do teste (χ², p-valor, MAD) ficam nos campos do struct, não na
# tabela — repeti-los em toda linha seria ruído.

Tables.istable(::Type{<:BenfordResult}) = true
Tables.rowaccess(::Type{<:BenfordResult}) = true
function Tables.rows(r::BenfordResult)
    n = r.n
    return [(digit = r.labels[i],
             count = r.counts[i],
             observed = r.observed[i],
             expected = r.expected[i],
             excess = r.observed[i] - r.expected[i],
             pearson_residual = (r.counts[i] - n * r.expected[i]) /
                                sqrt(n * r.expected[i]))
            for i in eachindex(r.labels)]
end

Tables.istable(::Type{<:LastDigitResult}) = true
Tables.rowaccess(::Type{<:LastDigitResult}) = true
function Tables.rows(r::LastDigitResult)
    e = 0.1
    return [(digit = d,
             count = r.counts[d + 1],
             observed = r.observed[d + 1],
             expected = e,
             excess = r.observed[d + 1] - e,
             pearson_residual = (r.counts[d + 1] - r.n * e) / sqrt(r.n * e))
            for d in 0:9]
end

Tables.istable(::Type{<:RozenasResult}) = true
Tables.rowaccess(::Type{<:RozenasResult}) = true
function Tables.rows(r::RozenasResult)
    return [(fraction = r.fractions[j],
             numerator = numerator(r.fractions[j]),
             denominator = denominator(r.fractions[j]),
             observed = r.observed[j],
             null_mean = r.null_mean[j],
             null_sd = r.null_sd[j],
             zscore = r.null_sd[j] > 0 ?
                      (r.observed[j] - r.null_mean[j]) / r.null_sd[j] :
                      (r.observed[j] > r.null_mean[j] ? Inf : 0.0),
             qvalue = r.qvalues[j])
            for j in eachindex(r.fractions)]
end

# `Tables.schema` explícito: sem isso o consumidor precisa inferir linha a linha.
Tables.schema(r::BenfordResult) = Tables.Schema(
    (:digit, :count, :observed, :expected, :excess, :pearson_residual),
    (Int, Int, Float64, Float64, Float64, Float64))
Tables.schema(r::LastDigitResult) = Tables.Schema(
    (:digit, :count, :observed, :expected, :excess, :pearson_residual),
    (Int, Int, Float64, Float64, Float64, Float64))
Tables.schema(r::RozenasResult) = Tables.Schema(
    (:fraction, :numerator, :denominator, :observed,
     :null_mean, :null_sd, :zscore, :qvalue),
    (Rational{Int}, Int, Int, Int, Float64, Float64, Float64, Float64))
