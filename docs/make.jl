push!(LOAD_PATH,"../src/")
using Documenter, SmoothedSpectralAbscissa
const SSA=SmoothedSpectralAbscissa

makedocs(;
    modules=[SmoothedSpectralAbscissa],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/dylanfesta/SmoothedSpectralAbscissa.jl/blob/{commit}{path}#L{line}",
    sitename="SmoothedSpectralAbscissa.jl",
    authors="Dylan Festa",
)

deploydocs(;
    repo="github.com/dylanfesta/SmoothedSpectralAbscissa.jl",
)
