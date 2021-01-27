using Documenter, DPOMPs

makedocs(sitename="DPOMPs.jl docs", pages = ["index.md", "models.md", "examples.md", "manual.md"])

## nb. called by GitHub Actions wf
# - local version deploys to build dir
deploydocs(
    repo = "github.com/mjb3/DPOMPs.jl.git",
)
