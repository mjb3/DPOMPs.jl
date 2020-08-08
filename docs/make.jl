using Documenter, DPOMPs

makedocs(sitename="DPOMPs.jl docs", pages = ["index.md", "models.md", "examples.md", "manual.md"])

deploydocs(
    repo = "github.com/mjb3/DPOMPs.jl.git",
)
