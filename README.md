# ElectionForensics.jl

[![CI](https://github.com/dantebertuzzi/ElectionForensics.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/dantebertuzzi/ElectionForensics.jl/actions/workflows/CI.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Election forensics in Julia: the coarse-fraction test of **Rozenas (2017)**,
**Benford**'s law (1st and 2nd digits), **last/penultimate digit** tests
(Beber & Scacco 2012), and the **election fingerprint** of Klimek et al.
(2012), with an integrated terminal report.

What sets this package apart is not the list of tests but its honesty about
when they hold: several classical nulls from the literature are **not** valid
for electoral counts, and `calibration_check` measures the type I error rate
each test actually delivers on *your* data before you conclude anything.

## Installation

```julia
] add https://github.com/dantebertuzzi/ElectionForensics.jl
```

## Quick start

```julia
using ElectionForensics

# votes[i]  = votes for the candidate in polling station i
# totals[i] = valid votes in polling station i

res = forensics_report(votes, totals)   # full battery, colored output

res.rozenas.pvalue
res.benford.conformity
```

Individual tests:

```julia
benford_test(votes; digit = 2)                    # 2BL, resampled null
last_digit_test(votes)                            # Beber & Scacco
last_digit_test(votes; position = :penultimate)   # ⚠ fragile uniform null
rozenas_test(votes, totals; boundary = :logit, B = 999)
rozenas_test(votes, totals; null = :betabinomial) # the `spikes` null
coarse_fractions(10)                              # 1/2, 1/3, 2/3, 1/4, ...
calibration_check(votes, totals)                  # does the p-value mean anything here?
```

Counts coming from CSV (`Float64`, `missing`) are coerced explicitly:
fractional values and `NaN` raise an error, and `missing` requires
`skipmissing = true`.

Five of the six result types implement Tables.jl — one row per digit, per
fraction, per test, or per occupied histogram cell:

```julia
using Tables, DataFrames
DataFrame(benford_test(votes))          # digit, count, observed, expected, excess, ...
DataFrame(rozenas_test(votes, totals))  # fraction, observed, null_mean, zscore, qvalue
```

## Methods

**Rozenas (2017).** Counts polling stations whose vote share `votes/totals` is
*exactly* an irreducible fraction k/d with d ≤ `max_denom` (compared in
rational arithmetic, with no floating-point error). The null distribution comes
from a parametric bootstrap: shares perturbed by a Gaussian kernel and counts
resampled from `Binomial(totals[i], p̃ᵢ)`. One-sided p-value with the
`(1 + #{T* ≥ T}) / (B + 1)` correction.

The jitter is applied on the **logit** scale (`boundary = :logit`, the
default): near 0 and 1 the step on the probability scale shrinks by itself, so
the null does not push mass from near-unanimous stations toward the middle,
where coarse fractions are dense. With jitter on the probability scale
(`:reflect`), under U-shaped shares the type I error fell to 0.000 and power at
ε = 2 % was 0.48 against 0.67.

Two nulls are available. `null = :kernel` (the default) perturbs the observed
shares. `null = :betabinomial` implements the model from Rozenas's own `spikes`
package: a mixture of Beta-Binomials fitted by EM and resampling from the
`Beta(y+α, n−y+β)` posterior, which deconvolves the binomial noise.

In simulation the two are equivalent. On real data the parametric one can
misfit a narrow band of shares and manufacture excess there — on the Russian
2012 data it underestimates the left tail by 3.8× and produces ten falsely
significant fractions. Every fitted mixture carries a `misfit` diagnostic, and
the package warns when it exceeds 2. Per-fraction excesses are exploratory: use
the `qvalues` (Benjamini–Hochberg), not the raw z-scores.

**Benford.** χ² goodness-of-fit plus MAD with Nigrini's (2012) conformity
thresholds. The default is the **2nd digit** (2BL, Mebane 2008) with
`null = :resampled`. Clean electoral counts do **not** follow Benford's law
when polling stations are homogeneous in size: against the classical null the
test rejects clean simulated elections in 72 % of cases (1BL: 100 %). The
resampled null — Gaussian jitter on `log10(count)` — brings the rate back to
0.05 and answers the useful question: "is the digit anomalous *given* the
empirical shape of the counts?". Use `null = :benford` as a descriptive
statistic only.

**Fingerprint (Klimek et al. 2012).** A 2-D histogram of polling stations by
turnout (x axis) and winner vote share (y axis). A single compact cluster is
the clean pattern; a smear toward the top-right corner indicates incremental
manipulation, and a second cluster at (100 %, 100 %) indicates extreme
manipulation. It also reports the skewness and kurtosis of the logarithmic vote
rate `ν = log((N−W)/W)`, which sit near (0, 3) in clean elections.

It requires one input more than the other tests — the **electorate** per
polling station:

```julia
election_fingerprint(winner_votes, valid_votes, electorate)
```

It is a visual diagnostic, not a test: it produces no p-value. The paper's
parametric fraud model `(fᵢ, f_e)` is not implemented — its specification lives
in the Supporting Information, outside the preprint.

**Last digit.** χ² against a uniform on 0:9, excluding small counts
(`min_value = 10`). Also reports `freq{0,5}` (expected ≈ 0.20; fabricated
numbers overuse 0 and 5). Calibrated for stations with ≥ 60 voters. The
**penultimate** digit requires the density to be smooth on a scale of 100 and
is therefore opt-in: with stations of 150–900 voters its uniform null rejects
clean elections in 50 % of cases (100 % with stations of 100–200). The
`null_valid` field flags when the condition is not met.

## Measure the calibration before concluding

The validity of each test's null depends on the shape of your count
distribution, and no simple statistic predicts when it fails. Rather than
promise, the package measures:

```julia
calibration_check(votes, totals)   # simulates clean elections with YOUR totals
```

It reports the type I error rate each test actually delivers on that data. The
verdict `:anticonservador` means the test rejects clean elections above α — do
not conclude fraud from it. `forensics_report` runs this check by default and
prints it before any p-value.

> **Note on language.** The public API, docstrings, and returned symbols
> (`:calibrado`, `:conservador`, `:anticonservador`) are in Portuguese, as is
> the terminal output. Only this README is in English.

## Caveats

No single test proves fraud. False positives arise from small polling stations,
administrative rounding, and aggregation; read the tests together and in
institutional context. Jittering around the observed shares makes the Rozenas
test slightly conservative under massive fraud — a rejection is therefore
strong evidence.

**What these tests do not detect.** None of them has power against
proportional ballot stuffing: in simulation, adding 15 % of votes to 20 % of
polling stations leaves all four tests at α. Rozenas detects round percentage
targets; the last-digit test detects decimal rounding. Fraud that leaves no
signature in digits or in exact fractions goes unnoticed.

The per-fraction excesses reported by the Rozenas test are **exploratory**:
under H₀ the largest z among the 31 fractions has median 2.4 and exceeds 2 in
68 % of clean elections. Use the `qvalues` (Benjamini–Hochberg), not the raw
z-scores.

## Validation

`validation/` contains the calibration and power scripts backing the claims
above. Run them with `julia -t auto validation/<script>.jl`.

## References

- Rozenas, A. (2017). Detecting Election Fraud from Irregularities in
  Vote-Share Distributions. *Political Analysis* 25(1), 41–56.
- Beber, B.; Scacco, A. (2012). What the Numbers Say: A Digit-Based Test
  for Election Fraud. *Political Analysis* 20(2), 211–234.
- Mebane, W. (2008). Election Forensics: The Second-Digit Benford's Law
  Test and Recent American Presidential Elections.
- Nigrini, M. (2012). *Benford's Law*. Wiley.
- Deckert, J.; Myagkov, M.; Ordeshook, P. (2011). Benford's Law and the
  Detection of Election Fraud. *Political Analysis* 19(3), 245–268.
- Klimek, P.; Yegorov, Y.; Hanel, R.; Thurner, S. (2012). Statistical
  Detection of Systematic Election Irregularities. *PNAS* 109(41),
  16469–16473.
- Benjamini, Y.; Hochberg, Y. (1995). Controlling the False Discovery Rate.
  *JRSS B* 57(1), 289–300.
- Silverman, B. (1986). *Density Estimation for Statistics and Data
  Analysis*. Chapman & Hall.

## License

MIT — see [LICENSE](LICENSE).
