module FWFTables

import Base.iterate
import Tables
import Parsers
import Formatting
import FixedSizeStrings

export readbla, Varspec, FWFTable, makefmt, File, write, stringtoint, CFILE

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
names(f::FWFTable) = [Symbol(var.name) for var in f.specs]
Tables.columnnames(f::FWFTable) = [Symbol(var.name) for var in f.specs]
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


abstract type AbstractNchar{N} <: AbstractString end
struct Nchar{N} <: AbstractNchar{N}
  value::String
end
Base.length(x::Nchar{N}) where N = N
Core.String(c::Nchar{N}) where N = c.value
Base.display(x::Nchar{N}) where N = display(x.value)
Base.show(io, x::Nchar{N}) where N = Base.show(io, x.value)


# Nchar is vervangen door FixedSizeString uit FixedSizeStrings
struct CharVector{N, L} <: AbstractVector{FixedSizeStrings.FixedSizeString{N}}
    buffer::Vector{UInt8}
    offset::Int64
    recordlength::Int64
end

Base.IndexStyle(cv::CharVector{N, L}) where {N, L} = IndexLinear()
Base.size(cv::CharVector{N, L}) where {N, L} = (L, 1)
function Base.getindex(cv::CharVector{N, L}, i::Int) where {N, L}
    startpos = (i-1) * cv.recordlength + cv.offset
    endpos = (i-1) * cv.recordlength + cv.offset + N - 1
    s = String(cv.buffer[startpos:endpos])
    return s
end

function Base.setindex!(cv::CharVector{N, L}, v::String, i::Int) where {N, L}
    startpos = (i-1) * cv.recordlength + cv.offset
    endpos = (i-1) * cv.recordlength + cv.offset + N - 1
    cv.buffer[startpos:endpos] = UInt8.(collect(v))
end

function Base.similar(cv::CharVector{N, L}) where {N, L}
    bufferlength = L * N
    return CharVector{N, L}(Vector{UInt8}(undef, bufferlength), 1, N)
end

function Base.copy(cv::CharVector{N, L}) where {N, L}
    bufferlength = L * N
    buffer = Vector{UInt8}(undef, bufferlength)
    offset = cv.offset
    recordlength = cv.recordlength
    for i=0:(L-1)
        buffer[i*N + 1:i*N + N] = cv.buffer[i*recordlength + offset: i*recordlength + offset + N - 1]
    end
    return CharVector{N, L}(buffer, 1, N)
end 

Tables.allocatecolumn(::Type{FixedSizeStrings.FixedSizeString{N}}, L) where N = CharVector{N, L}(Vector{UInt8}(undef, L*N), 1, N)
function Base.copyto!(dest::CharVector{N, L}, d_o::Integer, src::CharVector{N, M}) where {N, L, M}
    for i=0:(M-1)
        dest.buffer[(d_o-1+i)*dest.recordlength+dest.offset:(d_o-1+i)*dest.recordlength + dest.offset + N - 1] = src.buffer[i*src.recordlength + src.offset:i*src.recordlength + src.offset + N - 1]
    end
end

struct CFWFTable <: Tables.AbstractColumns
    specs::Vector{Varspec}
    columns::Dict{Symbol, AbstractVector}
end

Tables.istable(::CFWFTable) = true
Tables.columnaccess(::CFWFTable) = true

specs(t::CFWFTable) = getfield(t, :specs)
cols(t::CFWFTable) = getfield(t, :columns)

Tables.getcolumn(t::CFWFTable, i::Int) = cols(t)[Symbol(specs(t)[i].name)]
Tables.getcolumn(t::CFWFTable, nm::Symbol) = cols(t)[nm]
Tables.columnnames(t::CFWFTable) = [Symbol(spec.name) for spec in specs(t)]
names(t::CFWFTable) = [Symbol(var.name) for var in specs(t)]
types(t::CFWFTable) = [var.datatype for var in specs(t)]
Tables.schema(t::CFWFTable) = Tables.schema(names(t), types(t))
Base.length(t::CFWFTable) = length(cols(t)[Symbol(specs(t)[1].name)])


function CFile(filename::String, specs::Vector{Varspec})
    recordlength = maximum(spec.slice.stop for spec in specs)
    specselectie = [x for x in specs if x.datatype !== Nothing]
    buffer = open(filename) do io
        bom = read(io, 3) == [0xef, 0xbb, 0xbf]
        if !bom
            seekstart(io)
        end
        read(io)
    end
    # TODO geen records apart opvangen
    while buffer[recordlength+1] in [0x0a, 0x0d]
        recordlength = recordlength + 1
    end
    nrow = length(buffer) ÷ recordlength
    columns = Dict{Symbol, AbstractVector}()
    for spec in specselectie
        if spec.datatype == String
            length = spec.slice.stop - spec.slice.start + 1
            column = CharVector{length, nrow}(buffer, spec.slice.start, recordlength) 
        elseif spec.datatype == Union{Missing, Int64}
            column = Vector{Union{Missing, Int64}}(missing, nrow)
            for i in 1:nrow
                try
                    setindex!(column, bytestoint(buffer[(i-1)*recordlength + spec.slice.start:(i-1)*recordlength + spec.slice.stop]), i)
                catch e
                end
            end
        elseif spec.datatype == Float64
            column = Vector{Float64}(undef, nrow)
            for i in 1:nrow
                try
                    setindex!(column, bytestofloat(buffer[(i-1)*recordlength + spec.slice.start:(i-1)*recordlength + spec.slice.stop]), i)
                catch e
                    setindex!(column, NaN64, i)
                end
            end
        end
        columns[Symbol(spec.name)] = column
    end
    return CFWFTable(specselectie, columns)
end

function bytestoint(b)
    res = 0
    for c in b
        res = res * 10 + c - 48
    end
    return res
end

function bytestofloat(b, dec=0x2e)
    res = 0
    pre = true
    factor = 0.1
    for c in b
        if c == dec
            pre = false
            continue
        end
        if pre
            res = res * 10 + c - 48
        else
            res = res + factor * (c-48)
            factor = factor / 10
        end
    end
    return res
end


end # module

# vim: ts=4:sw=4:textwidth=89:fileencoding=utf-8
