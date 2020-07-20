using Test
using FWFTables
using FixedSizeStrings

bla = readbla("../testdata/testbla.bla")

@test length(bla) == 3
@test bla[1] == Varspec("var1", Union{Missing, Int64}, 1, 9, 0)
@test bla[2] == Varspec("var2", FixedSizeStrings.FixedSizeString, 10, 8, 0)
@test bla[3] == Varspec("var3", Float64, 18, 8, 2)

