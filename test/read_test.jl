using Test
using FWFTables
using InlineStrings


specs = [
         Varspec("var1", InlineStrings.String3, String, identity, 1, 2),
         Varspec("var2", InlineStrings.String3, String, identity, 3, 2),
        ]
tmpname = tempname()
Base.write(tmpname, "0102\n0304\n0506\n")
tb = FWFTables.File(tmpname, specs)

@test length(tb) == 3
@test tb.var1[1] == "01"
@test tb.var1[2] == "03"
@test tb.var1[3] == "05"
@test tb.var2[1] == "02"
@test tb.var2[2] == "04"
@test tb.var2[3] == "06"
                    
# no newline at end
tmpname = tempname()
Base.write(tmpname, "0102\n0304\n0506")
tb = FWFTables.File(tmpname, specs)
@test length(tb) == 3

# empty input
tmpname = tempname()
Base.write(tmpname, "")
tb = FWFTables.File(tmpname, specs)
@test length(tb) == 0

# empty input
tmpname = tempname()
Base.write(tmpname, "\n")
tb = FWFTables.File(tmpname, specs)
@test length(tb) == 0

