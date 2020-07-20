# julia -i -q --color=yes --project revise.jl example

using Revise, Jive
using FWFTables

trigger = function (path)
    printstyled("changed ", color=:cyan)
    println(path)
    revise()
    runtests(@__DIR__, skip=["revise.jl"])
end

watch(trigger, @__DIR__, sources=[pathof(FWFTables)])
trigger("")

Base.JLOptions().isinteractive==0 && wait()

