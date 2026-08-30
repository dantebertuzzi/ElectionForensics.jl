using Documenter
using ElectionForensics

DocMeta.setdocmeta!(ElectionForensics, :DocTestSetup,
                    :(using ElectionForensics); recursive = true)

makedocs(
    sitename = "ElectionForensics.jl",
    modules  = [ElectionForensics],
    authors  = "Dante Bertuzzi",
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical  = "https://dantebertuzzi.github.io/ElectionForensics.jl",
    ),
    pages = [
        "Início"     => "index.md",
        "Calibração" => "calibracao.md",
        "Referência" => "api.md",
    ],
    checkdocs = :exports,
    doctest   = true,
)

# comente `deploydocs` ao rodar localmente
deploydocs(
    repo = "github.com/dantebertuzzi/ElectionForensics.jl",
    devbranch = "main",
)
