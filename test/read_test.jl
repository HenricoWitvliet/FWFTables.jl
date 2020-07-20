using Test
using FWFTables
using FixedSizeStrings


specs = [
         Varspec("var1", FixedSizeString, 1, 2, 0),
         Varspec("var2", FixedSizeString, 3, 2, 0),
        ]
tb = FWFTables.File(IOBuffer("0102\n0304\n0506\n"), specs)

@test length(tb) == 3
@test tb.var1[1] == FixedSizeString{2}("01")
@test tb.var1[2] == FixedSizeString{2}("03")
@test tb.var1[3] == FixedSizeString{2}("05")
@test tb.var2[1] == FixedSizeString{2}("02")
@test tb.var2[2] == FixedSizeString{2}("04")
@test tb.var2[3] == FixedSizeString{2}("06")
                    
# no newline at end
tb = FWFTables.File(IOBuffer("0102\n0304\n0506"), specs)
@test length(tb) == 3

# empty input
tb = FWFTables.File(IOBuffer(""), specs)
@test length(tb) == 0

# empty input
tb = FWFTables.File(IOBuffer("\n"), specs)
@test length(tb) == 0

