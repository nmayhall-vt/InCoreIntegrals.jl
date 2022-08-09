using InCoreIntegrals
using Documenter

DocMeta.setdocmeta!(InCoreIntegrals, :DocTestSetup, :(using InCoreIntegrals); recursive=true)

makedocs(;
    modules=[InCoreIntegrals],
    authors="Nick Mayhall <nmayhall@vt.edu> and contributors",
    repo="https://github.com/nmayhall-vt/InCoreIntegrals.jl/blob/{commit}{path}#{line}",
    sitename="InCoreIntegrals.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://nmayhall-vt.github.io/InCoreIntegrals.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/nmayhall-vt/InCoreIntegrals.jl",
    devbranch="main",
)
