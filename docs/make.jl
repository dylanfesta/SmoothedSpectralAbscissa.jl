push!(LOAD_PATH,abspath(@__DIR__,"..","..","src"))
using Documenter
using SmoothedSpectralAbscissa ; const SSA=SmoothedSpectralAbscissa
using Plots,NamedColors
using LinearAlgebra
using Random

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
    devbranch="main",
    devurl="dev",
    versions = ["stable" => "v^", "v#.#", devurl => devurl],
)
