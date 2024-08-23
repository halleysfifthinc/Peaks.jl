using Documenter, DocumenterInterLinks
using Peaks

ENV["GKSwstype"] = "100"

DocMeta.setdocmeta!(Peaks, :DocTestSetup, :(using Peaks); recursive=true)
makedocs(
    sitename = "Peaks",
    modules = [Peaks],
    checkdocs = :exports,
    format = Documenter.HTML(;
        prettyurls=true,
        canonical="https://halleysfifthinc.github.io/Peaks.jl",
    ),
    authors="Allen Hill <allenofthehills@gmail.com> and contributors",
    repo=Remotes.GitHub("halleysfifthinc", "Peaks.jl"),
    pages=[
        "Home" => "index.md",
        "How-to" => "how-to.md",
        "Benchmarks" => "benchmarks.md",
        "Reference" => [
            "Glossary" => "glossary.md",
            "API" => "reference.md",
            ],
    ],
    plugins = [
        InterLinks(
            "Base" => "https://docs.julialang.org/en/v1/objects.inv"
        )
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/halleysfifthinc/Peaks.jl",
    push_preview=true,
)
