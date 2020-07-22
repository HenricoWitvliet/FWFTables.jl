FWFTables.jl
============

Deze bibliotheek bevat functies om eenvoudig fixed width bestanden te lezen en schrijven.

```julia
using DataFrames, FWFTables

df = DataFrame(FWFTables.File("databestand.asc", "spec.bla")
```
