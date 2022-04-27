FWFTables.jl
============

![GitHub CI](https://github.com/HenricoWitvliet/FWFTables.jl/actions/workflows/ci.yml/badge.svg)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://HenricoWitvliet.github.io/FWFTables.jl/dev)

This package contains functions for working with fixed width ascii files.

The File-function is used to read data. The specifications must be given as
Blaise-bla files.
The write-function will write a tables-object to a fixed width ascii file. The
specs must also be given as a bla-file.

If the file starts with a bom for utf-8, it is skipped. The File-function
tries to guess the used line-ending.


```julia
using DataFrames, FWFTables

df = DataFrame(FWFTables.File("databestand.asc", "spec.bla"))
```
