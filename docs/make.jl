using Documenter, DPOMPs

makedocs(sitename="DPOMPs.jl docs", doctest = false)

deploydocs(
    repo = "github.com/mjb3/DPOMPs.jl.git",
)
