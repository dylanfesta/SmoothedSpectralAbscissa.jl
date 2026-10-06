using Documenter
using Literate
using SmoothedSpectralAbscissa

const SSA = SmoothedSpectralAbscissa
DocMeta.setdocmeta!(SSA, :DocTestSetup, :(using SmoothedSpectralAbscissa); recursive=true)

for name in ("01_show_ssa", "02_dynamics", "03_EI")
    Literate.markdown(
        joinpath(@__DIR__, "..", "examples", name * ".jl"),
        joinpath(@__DIR__, "src", "generated");
        flavor=Literate.DocumenterFlavor(),
    )
end

makedocs(;
    modules=[SSA],
    authors="Dylan Festa",
    sitename="SmoothedSpectralAbscissa.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://dylanfesta.github.io/SmoothedSpectralAbscissa.jl",
        repolink="https://github.com/dylanfesta/SmoothedSpectralAbscissa.jl",
        edit_link="main",
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => [
            "Comparison of SA and SSA" => "generated/01_show_ssa.md",
            "Stability-optimized linear systems" => "generated/02_dynamics.md",
            "Excitatory/inhibitory recurrent network" => "generated/03_EI.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/dylanfesta/SmoothedSpectralAbscissa.jl",
    devbranch="main",
    devurl="dev",
    versions=["stable" => "v^", "v#.#", "dev" => "dev"],
)
