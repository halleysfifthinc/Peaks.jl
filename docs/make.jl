using Documenter
using Peaks

DocMeta.setdocmeta!(Peaks, :DocTestSetup, :(using Peaks); recursive=true)
makedocs(
    sitename = "Peaks",
    modules = [Peaks],
    checkdocs = :exports,
    format = Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://halleysfifthinc.github.io/Peaks.jl",
    ),
    authors="Allen Hill <allenofthehills@gmail.com> and contributors",
    repo=Remotes.GitHub("halleysfifthinc", "Peaks.jl")
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = Remotes.GitHub("halleysfifthinc", "Peaks.jl"),
    push_preview=true,
)
