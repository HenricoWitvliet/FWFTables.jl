module FWFTables

import Base.iterate
import Tables
import Parsers
import Formatting

export readbla, Varspec, FWFTable, makefmt, File, write, stringtoint

"""
    Blavar

Struct to hold the specification for a single variable from a Blaise specification.

# Examples
```julia-repl
julia> blavar = Blavar("var1", 1:9, String, 9, 0)
```
"""
struct Blavar
    name::String
    slice::UnitRange{Int64}
    datatype::Type
    length::Int64
    decimals::Int64
end

struct Varspec
    name::String
    slice::UnitRange{Int64}
    datatype::Type
    stringtodata::Function
    datatostring::Function
end


varregex = r"(?i)\s*((?P<name>\w+)\s*(?P<tekst>\".*\"|)\s*:|)\s*(?P<array>ARRAY\[(?P<astart>\d+)\.\.(?P<aeind>\d+)\]\s*OF\s*|)(?P<type>DUMMY|STRING|REAL|INTEGER)\s*\[\s*(?P<length>\d+)(\s*,\s*(?P<decimals>\d+)|)\s*\].*"

convdatatype = Dict(
    "string" => String,
    "integer" => Union{Missing,Int64},
    "dummy" => Nothing,
    "real" => Float64,
)

Parsers.tryparse(::Type{Int}, ::Nothing) = 0
Parsers.tryparse(::Type{String}, s::String) = s

stringtodummy(s::String) = nothing
stringtostring(s::String)::String = s
function stringtoint(s::String)::Union{Missing,Int64}
    val = try
        Parsers.parse(Int64, s)
    catch e
        missing
    end
    return val
end
function stringtofloat(s::String)::Float64
    val = try
        Parsers.parse(Float64, s)
    catch e
        NaN
    end
    return val
end

convdatafunctie = Dict(
    "string" => stringtostring,
    "integer" => stringtoint,
    "dummy" => stringtodummy,
    "real" => stringtofloat,
)


"""
    readbla(filename)

Read a Blaise specification file for a fixed width ascii file. Returns a vector
of 'Varspec' objects.
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
            decimals = Parsers.tryparse(Int64, var[:decimals])
            datatype = convdatatype[lowercase(var[:type])]
            slice = startpos:(startpos+len-1)
            blavar = Blavar(name, slice, datatype, len, decimals)
            dataconv = convdatafunctie[lowercase(var[:type])]
            datafmt = makefmt(blavar)
            if var[:array] == ""
                push!(spec, Varspec(name, slice, datatype, dataconv, datafmt))
                startpos += len
            else
                # vector of vars
                for varnr = parse(Int64, var[:astart]):parse(Int64, var[a:eind])
                    tmpname = name * "_" * string(varnr)
                    slice = startpos:(startpos+len-1)
                    blavar = Blavar(tmpname, slice, datatype, len, decimals)
                    push!(spec, Varspec(tmpname, slice, datatype, dataconv, datafmt))
                    startpos += len
                end
            end
        end
    end
    return spec
end

function makefmt(blavar::Blavar)
    if blavar.datatype == String
        fmt = "{:0" * string(blavar.length) * "s}"
    elseif blavar.datatype == Int64
        fmt = "{:0" * string(blavar.length) * "d}"
    elseif blavar.datatype == Float64
        fmt = "{:0" * string(blavar.length) * "." * string(blavar.decimals) * "f}"
    elseif blavar.datatype == Nothing
        fmt = repeat(" ", string(blavar.length))
    end
    return x -> fmt
end

function makefmt(specs::Vector{Varspec})
    string([spec.datatostring for spec in specs]...)
end

struct FWFTable
    handle::IOStream
    specs::Vector{Varspec}
    names::Dict{Symbol,Varspec}
    numberofrecords::Int
end

struct FWFTableRow <: Tables.AbstractRow
    row::Int
    rawrow::String
    source::FWFTable
end


function nextline(f, st)
    rawrow = readline(f.handle)
    if rawrow == ""
        return nothing
    else
        return (FWFTableRow(st + 1, rawrow, f), st + 1)
    end
end


Tables.istable(::FWFTable) = true
names(f::FWFTable) = [var.name for var in f.specs]
Tables.columnnames(f::FWFTable) = [var.name for var in f.specs]
types(f::FWFTable) = [var.datatype for var in f.specs]
Tables.schema(f::FWFTable) = Tables.Schema(names(f), types(f))


Tables.rowaccess(::FWFTable) = true
Tables.rows(f::FWFTable) = f
Base.eltype(f::FWFTable) = FWFTableRow
Base.length(f::FWFTable) = getfield(f, :numberofrecords)
Base.iterate(f::FWFTable, st = 0) = nextline(f, st)

Tables.getcolumn(r::FWFTableRow, s::String) = Tables.getcolumn(r, Symbol(s))

Tables.columnnames(r::FWFTableRow) = names(getfield(r, :source))

function Tables.getcolumn(r::FWFTableRow, col::Int)
    var = getfield(getfield(r, :source), :specs)[col]
    value = var.stringtodata(getfield(r, :rawrow)[var.slice])
    return value
end

function Tables.getcolumn(r::FWFTableRow, nm::Symbol)
    var = getfield(getfield(r, :source), :names)[nm]
    value = var.stringtodata(getfield(r, :rawrow)[var.slice])
    return value
end

"""
    File(filename, blafilename)

Read a fixed width file, using the specs from blafile. Returns an object
implementing the Tables interface. 

# Examples
```julia-repl
julia> using DataFrames, FWFTables
julia> df = DataFrame(FWFTables.File("data.asc", "spec.bla")
```
"""
function File(filename::String, blafilename::String)
    specs = readbla(blafilename)
    File(filename, specs)
end

"""
    File(filename, specs)

Use a vector of Varspec-definitions instead of the bla-file

# Examples
```julia-repl
julia> using DataFrames, FWFTables
julia> specs = readbla("spec.bla")
julia> df = DataFrame(FWFTables.File("data.asc", specs)
```
"""
function File(filename::String, specs::Vector{Varspec})
    specselectie = [x for x in specs if x.datatype !== Nothing]
    d = Dict([Symbol(elt.name) => elt for elt in specselectie])
    handle = open(filename)
    numberofrecords = countlines(handle)
    seekstart(handle)
    bom = read(handle, 3) == [0xef, 0xbb, 0xbf]
    if !bom
        seekstart(handle)
    end
    return FWFTable(handle, specselectie, d, numberofrecords)
end

function write(filename::String, blafilename::String, table)
    specs::Vector{Varspec} = readbla(blafilename)
    write(filename, specs, table)
end

function write(filename::String, specs::Vector{Varspec}, table)
    fe = Formatting.FormatExpr(makefmt(specs))
    open(filename, "w") do io
        for row in Tables.rows(table)
            Formatting.printfmtln(io, fe, row...)
        end
    end
end

end # module

# vim: ts=2:sw=2:textwidth=89
