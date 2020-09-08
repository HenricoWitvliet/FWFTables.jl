var documenterSearchIndex = {"docs":
[{"location":"#FWFTables.jl","page":"FWFTables.jl","title":"FWFTables.jl","text":"","category":"section"},{"location":"","page":"FWFTables.jl","title":"FWFTables.jl","text":"Documentation for FWFTables.jl","category":"page"},{"location":"","page":"FWFTables.jl","title":"FWFTables.jl","text":"Modules = [FWFTables]\nOrder   = [:type, :function]","category":"page"},{"location":"#FWFTables.Varspec","page":"FWFTables.jl","title":"FWFTables.Varspec","text":"Varspec\n\nStruct to hold all relevant information for one variable.  Fields are\n\nname: name of the variable\ndatatype: Julia type of variable\nstartpos: starting position of variabele in record\nlength: length in bytes of variable\ndecimals: for float-type the number of decimals\n\nPredefined usable types are FixedSizeString{N} for string variables, Union{Missing,Int64} and IntXX for integer variables, FloatXX for real numbers and Nothing for variables that must be skipped. If you want to use another type for input, then you must define fw_convert to convert the byte-string to the given type.\n\nTo write to a fixed width file, generator functions for the used types must be available. Predefined function are defined for FixedSizeString{N}, Union{Missing,Int64}, Float64 and Nothing. For other functions the function makewrite must be defined.\n\n\n\n\n\n","category":"type"},{"location":"#FWFTables.File-Tuple{IO,Array{Varspec,1}}","page":"FWFTables.jl","title":"FWFTables.File","text":"File(io:IO, specs::Vector{Varspec})\n\nRead a fixed width file.\n\n\n\n\n\n","category":"method"},{"location":"#FWFTables.File-Tuple{String,Array{Varspec,1}}","page":"FWFTables.jl","title":"FWFTables.File","text":"File(filename::String, specs::Vector{Varspec})\n\nRead a fixed width file.\n\n\n\n\n\n","category":"method"},{"location":"#FWFTables.File-Tuple{String,String}","page":"FWFTables.jl","title":"FWFTables.File","text":"File(filename::String, blafilename::String)\n\nRead a fixed width file, using the specs from blafile. Returns an object implementing the Tables interface. \n\nIf the file starts with the 3 byte BOM for utf-8, it will be skipped. The line-ending is determined based on the first line. If the file uses a mixed line ending, the data will not be read correctly.\n\nExamples\n\njulia> using DataFrames, FWFTables\njulia> df = DataFrame(FWFTables.File(\"data.asc\", \"spec.bla\"))\n\n\n\n\n\n","category":"method"},{"location":"#FWFTables.readbla-Tuple{Any}","page":"FWFTables.jl","title":"FWFTables.readbla","text":"readbla(filename)\n\nRead a Blaise specification file for a fixed width ascii file. Returns a vector of 'Varspec' objects.\n\nA simple example of a Blaise specification file looks like\n\ndatamodel testbla\n Fields\n    var1        : integer[9]\n    var2        : string[8]\n    dummy[2]\n    var3        : Real[8, 2]\n    vars        : Array[2010..2050] of STRING[4]\nendmodel\n\n\n\n\n\n","category":"method"},{"location":"#FWFTables.write-Tuple{IO,Array{Varspec,1},Any}","page":"FWFTables.jl","title":"FWFTables.write","text":"write(io::IO, specs::Vector{Varspec}, table)\n\nWrite the table to a fixed width ascii file\n\n\n\n\n\n","category":"method"},{"location":"#FWFTables.write-Tuple{String,Array{Varspec,1},Any}","page":"FWFTables.jl","title":"FWFTables.write","text":"write(filename::String, specs::Vector{Varspec}, table)\n\nWrite the table to a fixed width ascii file\n\n\n\n\n\n","category":"method"},{"location":"#FWFTables.write-Tuple{String,String,Any}","page":"FWFTables.jl","title":"FWFTables.write","text":"write(filename::String, blafilename::String, table)\n\nWrite the table to a fixed width ascii file. The blafile contains the specifications of the file in the Blaise format.\n\nUsing the datatype, length and if relevant the decimals from the specification a conversion function is generated to write a single value using the required width. For String, Integer, Real and Dummy there are conversion-functions generated by makewrite.  For other datatypes, the required generator must be created. The signature is (T::Type, spec::Varspec), and the return value must be a function that accepts (io::IO, value::T) and writes the ascii-representation of the value to io.\n\n\n\n\n\n","category":"method"}]
}