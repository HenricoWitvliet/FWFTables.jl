push!(LOAD_PATH, "../")

using Jive
runtests(@__DIR__, skip=["revise.jl"])

