using Test
using FWFTables
using InlineStrings

bla = readbla("../testdata/testbla.bla")

@test length(bla) == 3
b1 = bla[1]
b2 = bla[2]
b3 = bla[3]

@test b1.name == "var1"
@test b1.storagetype == Union{Missing, Int64}
@test b1.startpos == 1
@test b1.length == 9
@test b2.name == "var2"
@test b2.storagetype == InlineStrings.String15
@test b2.startpos == 10
@test b2.length == 8
@test b3.name == "var3"
@test b3.storagetype == Float64
@test b3.startpos == 18
@test b3.length == 8

