push!(LOAD_PATH,"../src/")
using Documenter, SmoothedSpectralAbscissa
const SSA=SmoothedSpectralAbscissa

DocMeta.setdocmeta!(SmoothedSpectralAbscissa, :DocTestSetup, :(using SmoothedSpectralAbscissa); recursive=true)

makedocs(;
    modules=[SmoothedSpectralAbscissa],
    authors="Dylan Festa",
    repo="https://github.com/dylanfesta/SmoothedSpectralAbscissa.jl/blob/{commit}{path}#{line}",
    sitename="SmoothedSpectralAbscissa.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://dylanfesta.github.io/SmoothedSpectralAbscissa.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/dylanfesta/SmoothedSpectralAbscissa.jl",
    devbranch="master",
)
