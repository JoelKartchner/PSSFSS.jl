using PSSFSS
using Documenter
using DemoCards

#=
using Literate
cd("literate") do
  include("literate/compile.jl") 
end
=#

demopage, postprocess_cb, demo_assets = makedemos("PSS_&_FSS_Element_Gallery")
assets = String[]
isnothing(demo_assets) || (push!(assets, demo_assets))

makedocs(;
    clean=false,
    modules=[PSSFSS],
    authors="Peter Simon <psimon0420@gmail.com> and contributors",
    repo="https://github.com/simonp0420/PSSFSS.jl/blob/{commit}{path}#L{line}",
    sitename="PSSFSS.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://simonp0420.github.io/PSSFSS.jl/stable",
        assets=assets,
    ),
    pages=[
        "Home" => "index.md",
        "User Manual" => "manual.md",
        demopage,
        "Usage Examples" => "examples.md",
        "Function Reference" => "reference.md",
        "Index" => "function_index.md"
    ],
)
postprocess_cb() # For DemoCards

deploydocs(;
    repo="github.com/simonp0420/PSSFSS.jl",
    devbranch = "main"
)
