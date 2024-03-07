module FWFTables

import Tables
import Parsers
import Format
import InlineStrings
import Mmap

export readbla, Varspec, FWFTable, File, write

"""
    Varspec

Struct to hold all relevant information for one variable. 
Fields are

- name: name of the variable
- storagetype: Julia type of variable
- conversionfcie: Julia function to convert the raw string value to the Julia type, used for reading
- dumpfcie: Julia function to convert the Julia type to a fixed width string, used for writing
- startpos: starting position of variabele in record
- length: length in bytes of variable

Predefined usable types are `InlineStringXX`-types for string variables,
`Union{Missing,Int64}` and `IntXX` for integer variables, `FloatXX` for real numbers and
`Nothing` for variables that must be skipped. If you want to use another type for input,
then you must specify a function to convert the byte-string to the given type.

To write to a fixed width file, conversion functions for the used types must be available.
Predefined function are defined for `String`, `Union{Missing,Int64}`,
`Float64` and `Nothing`.
"""
struct Varspec
    name::String
    storagetype::Type
    conversionfcie::Union{DataType,Function}
    dumpfcie::Union{DataType,Function}
    startpos::Int64
    length::Int64
end

varregex = r"(?i)\s*((?P<name>\w+)\s*(?P<tekst>\".*\"|)\s*:|)\s*(?P<array>ARRAY\[(?P<astart>\d+)\.\.(?P<aeind>\d+)\]\s*OF\s*|)(?P<type>DUMMY|STRING|REAL|INTEGER)\s*\[\s*(?P<length>\d+)(\s*,\s*(?P<decimals>\d+)|)\s*\].*"

convdatatype = Dict(
    "string" => String,
    "integer" => Union{Missing,Int64},
    "dummy" => Nothing,
    "real" => Float64,
)

convint = value -> Parsers.tryparse(Int64, value)
convfloat = value -> something(Parsers.tryparse(Float64, value), NaN64)

convfcie = Dict(
    "string" => identity,
    "integer" => convint,
    "dummy" => identity,
    "real" => convfloat
)

Base.tryparse(::Type{Int}, ::Nothing) = 0


"""
    readbla(filename)

Read a Blaise specification file for a fixed width ascii file. Returns a vector
of 'Varspec' objects.

A simple example of a Blaise specification file looks like
```
datamodel testbla
 Fields
    var1        : integer[9]
    var2        : string[8]
    dummy[2]
    var3        : Real[8, 2]
    vars        : Array[2010..2050] of STRING[4]
endmodel
```

"""
function readbla(filename)
    spec = Vector{Varspec}()
    regels = open(filename) do io
        readlines(io)
    end
    parsevars = false
    startpos = 1
    for regel in regels
        if !parsevars
            if occursin(r"(?i)fields", regel)
                parsevars = true
            end
            continue
        end
        var = match(varregex, regel)
        if !isnothing(var)
            if !isnothing(var[:name])
                name = lowercase(var[:name])
            else
                name = "dummy"
            end
            len = parse(Int64, var[:length])
            decimals = Base.tryparse(Int64, var[:decimals])
            storagetype = convdatatype[lowercase(var[:type])]
            fcie = convfcie[lowercase(var[:type])]
            if storagetype == String
                storagetype = len > 255 ? String : InlineStrings.InlineStringType(len)
                fcie = storagetype
                fe = Format.FormatSpec(">" * string(len) * "s")
                dumpfcie = (io, value) -> Format.printfmt(io, fe, value)
            else
                dumpfcie = makewrite(storagetype, len, decimals)
            end
            if var[:array] == ""
                push!(spec, Varspec(name, storagetype, fcie, dumpfcie, startpos, len))
                startpos += len
            else
                # vector of vars
                for varnr = parse(Int64, var[:astart]):parse(Int64, var[:aeind])
                    tmpname = name * "_" * string(varnr)
                    push!(spec, Varspec(tmpname, storagetype, fcie, dumpfcie, startpos, len))
                    startpos += len
                end
            end
        end
    end
    return spec
end


struct FWFTable <: Tables.AbstractColumns
    specs::Vector{Varspec}
    columns::Dict{Symbol,AbstractVector}
end

Tables.istable(::FWFTable) = true
Tables.columnaccess(::FWFTable) = true

specs(t::FWFTable) = getfield(t, :specs)
cols(t::FWFTable) = getfield(t, :columns)

Tables.getcolumn(t::FWFTable, i::Int) = cols(t)[Symbol(specs(t)[i].name)]
Tables.getcolumn(t::FWFTable, nm::Symbol) = cols(t)[nm]
Tables.columnnames(t::FWFTable) = [Symbol(spec.name) for spec in specs(t)]
names(t::FWFTable) = [Symbol(var.name) for var in specs(t)]
types(t::FWFTable) = [var.storagetype for var in specs(t)]
Tables.schema(t::FWFTable) = Tables.Schema(names(t), types(t))
Base.length(t::FWFTable) = length(cols(t)[Symbol(specs(t)[1].name)])

function createcolumn(T, fcie, buffer, nrow, startpos, length, recordlength)
    column = Vector{T}(undef, nrow)
    pos = startpos
    endpos = startpos + length - 1
    for i = 1:nrow
        setindex!(column, fcie(view(buffer, pos:endpos)), i)
        pos = pos + recordlength
        endpos = endpos + recordlength
    end
    return column
end


"""
    File(filename::String, blafilename::String)

Read a fixed width file, using the specs from blafile. Returns an object
implementing the Tables interface. 

If the file starts with the 3 byte BOM for utf-8, it will be skipped.
The line-ending is determined based on the first line. If the file uses a mixed line
ending, the data will not be read correctly.

# Examples
```julia-repl
julia> using DataFrames, FWFTables
julia> df = DataFrame(FWFTables.File("data.asc", "spec.bla"))
```
"""
function File(filename::String, blafilename::String)
    specs = readbla(blafilename)
    File(filename, specs)
end

"""
    File(filename::String, specs::Vector{Varspec})

Read a fixed width file.
"""
function File(filename::String, specs::Vector{Varspec})
    open(filename) do io
        File(io, specs, 2^62)
    end
end


function File(filename::String, specs::Vector{Varspec}, startrow::Int64, nrow::Int64)
    recordlength = maximum(spec.startpos + spec.length - 1 for spec in specs)
    open(filename) do io
        pos = position(io)
        bom = read(io, 3) == [0xef, 0xbb, 0xbf]
        if !bom
            seek(io, pos)
        end
        buffer = Mmap.mmap(io)
        rawlength = recordlength
        while (length(buffer) > recordlength) && (buffer[recordlength+1] in [0x0a, 0x0d])
            recordlength = recordlength + 1
        end
        crlflength = recordlength - rawlength
        # calculate max number of rows based on record length
        # make a correction if last line doesn't contain eol
        if crlflength > 0 && !(buffer[end] in [0x0a, 0x0d])
            nrowmax = (length(buffer) + crlflength) ÷ recordlength
        else
            nrowmax = length(buffer) ÷ recordlength
        end
        # 
        if startrow <= nrowmax
            startpos = (bom ? 3 : 0) + (startrow - 1) * recordlength
            seek(io, startpos)
            File(io, specs, nrow)
        end
    end

end


"""
    File(io:IO, specs::Vector{Varspec})

Read a fixed width file.
"""
function File(io::IO, specs::Vector{Varspec})
    return(File(io, specs, 2^62))
end


function File(io::IO, specs::Vector{Varspec}, nrow::Int64)
    recordlength = maximum(spec.startpos + spec.length - 1 for spec in specs)
    specselectie = [x for x in specs if x.storagetype !== Nothing]
    pos = position(io)
    bom = read(io, 3) == [0xef, 0xbb, 0xbf]
    if !bom
        seek(io, pos)
    end
    buffer = Mmap.mmap(io)
    # read a line to determine line ending
    rawlength = recordlength
    while (length(buffer) > recordlength) && (buffer[recordlength+1] in [0x0a, 0x0d])
        recordlength = recordlength + 1
    end
    crlflength = recordlength - rawlength
    # calculate max number of rows based on record length
    # make a correction if last line doesn't contain eol
    if crlflength > 0 && !(buffer[end] in [0x0a, 0x0d])
        nrowmax = (length(buffer) + crlflength) ÷ recordlength
    else
        nrowmax = length(buffer) ÷ recordlength
    end
    nrow = min(nrow, nrowmax)
    columns = Dict{Symbol,AbstractVector}()
    Threads.@threads for spec in specselectie
        column = createcolumn(spec.storagetype, spec.conversionfcie, buffer, nrow, spec.startpos, spec.length, recordlength)
        columns[Symbol(spec.name)] = column
    end
    return FWFTable(specselectie, columns)
end


function makewrite(::Type{Union{Missing, Int64}}, length, decimals)
    fe = Format.FormatSpec(string("0>", length, "d"))
    spaces = repeat(" ", length) 
    return (io, val) -> ismissing(val) ? spaces : Format.printfmt(io, fe, val)
end
function makewrite(::Type{Float64}, length, decimals)
    fe = Format.FormatSpec(string("0>", length, ".", decimals, "f"))
    return (io, val) -> Format.printfmt(io, fe, val)
end
function makewrite(::Type{Nothing}, length, decimals)
    spaces = repeat(" ", length) 
    return io -> Base.write(io, spaces)
end


"""
    write(filename::String, blafilename::String, table)

Write the table to a fixed width ascii file. The blafile contains the specifications of
the file in the Blaise format.

Using the datatype, length and if relevant the decimals from the specification a conversion
function is generated to write a single value using the required width.
For String, Integer, Real and Dummy there are conversion-functions generated by
*makewrite*.  For other datatypes, the required generator must be provided.
"""
function write(filename::String, blafilename::String, table)
    specs::Vector{Varspec} = readbla(blafilename)
    write(filename, specs, table)
end

"""
    write(filename::String, specs::Vector{Varspec}, table)

Write the table to a fixed width ascii file
"""
function write(filename::String, specs::Vector{Varspec}, table)
    open(filename, "w") do io
        write(io, specs, table)
    end
end

"""
    write(io::IO, specs::Vector{Varspec}, table)

Write the table to a fixed width ascii file
"""
function write(io::IO, specs::Vector{Varspec}, table)
    writefcies = [Symbol(spec.name) => spec.dumpfcie for spec in specs]
    for row in Tables.rows(table)
        for (name, writefcie) in writefcies
            if name == :dummy
                writefcie(io)
            else
                writefcie(io, row[name])
            end
        end    
        println(io, "")
    end
end


end # module

# vim: ts=4:sw=4:textwidth=89:fileencoding=utf-8
