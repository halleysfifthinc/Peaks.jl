using Documenter
using Peaks

DocMeta.setdocmeta!(Peaks, :DocTestSetup, :(using Peaks); recursive=true)
makedocs(
    sitename = "Peaks",
    modules = [Peaks],
    format = Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://halleysfifthinc.github.io/Peaks.jl",
    ),
    authors="Allen Hill <allenofthehills@gmail.com> and contributors",
    repo="https://github.com/halleysfifthinc/Peaks.jl/blob/{commit}{path}#L{line}",
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/halleysfifthinc/Peaks.jl",
    push_preview=true,
)
