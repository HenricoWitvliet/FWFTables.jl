push!(LOAD_PATH,"../src/")
using Documenter
using FWFTables

makedocs(
    sitename = "FWFTables",
    format = Documenter.HTML(),
    modules = [FWFTables]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/HenricoWitvliet/FWFTables.jl"
)
