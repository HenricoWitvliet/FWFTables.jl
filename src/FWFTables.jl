module FWFTables

import Base.iterate
import Tables
import Parsers
import Formatting
import FixedSizeStrings.FixedSizeString

export readbla, Varspec, FWFTable, makefmt, File, write, stringtoint


struct Varspec
    name::String
    slice::UnitRange{Int64}
    datatype::Type
    startpos::Int64
    length::Int64
    decimals::Int64
end


varregex = r"(?i)\s*((?P<name>\w+)\s*(?P<tekst>\".*\"|)\s*:|)\s*(?P<array>ARRAY\[(?P<astart>\d+)\.\.(?P<aeind>\d+)\]\s*OF\s*|)(?P<type>DUMMY|STRING|REAL|INTEGER)\s*\[\s*(?P<length>\d+)(\s*,\s*(?P<decimals>\d+)|)\s*\].*"

convdatatype = Dict(
    "string" => FixedSizeString,
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
            if var[:array] == ""
                push!(spec, Varspec(name, slice, datatype, startpos, len, decimals))
                startpos += len
            else
                # vector of vars
                for varnr = parse(Int64, var[:astart]):parse(Int64, var[a:eind])
                    tmpname = name * "_" * string(varnr)
                    slice = startpos:(startpos+len-1)
                    push!(spec, Varspec(tmpname, slice, datatype, startpos, len, decimals))
                    startpos += len
                end
            end
        end
    end
    return spec
end

function makefmt(spec::Varspec)
    if spec.datatype == FixedSizeString
        fmt = "{:0" * string(spec.length) * "s}"
    elseif spec.datatype == Union{Missing, Int64}
        fmt = "{:0" * string(spec.length) * "d}"
    elseif spec.datatype == Float64
        fmt = "{:0" * string(spec.length) * "." * string(spec.decimals) * "f}"
    elseif spec.datatype == Nothing
        fmt = repeat(" ", string(spec.length))
    end
    return fmt
end

function makefmt(specs::Vector{Varspec})
    string([makefmt(spec) for spec in specs]...)
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
Base.length(x::Nchar{N}) where {N} = N
Core.String(c::Nchar{N}) where {N} = c.value
Base.display(x::Nchar{N}) where {N} = display(x.value)
Base.show(io, x::Nchar{N}) where {N} = Base.show(io, x.value)


# Nchar is vervangen door FixedSizeString uit FixedSizeStrings
struct CharVector{N,L} <: AbstractVector{FixedSizeString{N}}
    buffer::Vector{UInt8}
    offset::Int64
    recordlength::Int64
end

Base.IndexStyle(cv::CharVector{N,L}) where {N,L} = IndexLinear()
Base.size(cv::CharVector{N,L}) where {N,L} = (L, 1)

function Base.getindex(cv::CharVector{N,L}, i::Integer) where {N,L}
    startpos = (i - 1) * cv.recordlength + cv.offset
    endpos = (i - 1) * cv.recordlength + cv.offset + N - 1
    #s = String(cv.buffer[startpos:endpos])
    s = FixedSizeString{N}(view(cv.buffer, startpos:endpos))
    return s
end

Base.view(cv::CharVector{N, L}, i::Integer) where {N, L} = cv[i]
Base.getindex(cv::CharVector{N, L}, inds::AbstractUnitRange) where {N, L} = view(cv, collect(inds))
Base.getindex(cv::CharVector{N, L}, inds) where {N, L} = view(cv, inds)
Base.Broadcast.dotview(cv::CharVector{N, L}, inds::AbstractUnitRange) where {N, L} = Base.Broadcast.dotview(cv, collect(inds))

function Base.setindex!(cv::CharVector{N,L}, v, i::Integer) where {N,L}
    startpos = (i - 1) * cv.recordlength + cv.offset
    endpos = (i - 1) * cv.recordlength + cv.offset + N - 1
    cv.buffer[startpos:endpos] = UInt8.(collect(v))
end

function Base.setindex!(cv::CharVector{N,L}, v, inds::AbstractUnitRange) where {N,L}
    for i in inds
        startpos = (i - 1) * cv.recordlength + cv.offset
        endpos = (i - 1) * cv.recordlength + cv.offset + N - 1
        cv.buffer[startpos:endpos] = UInt8.(collect(v))
    end
end


function Base.similar(cv::CharVector{N,L}) where {N,L}
    bufferlength = L * N
    return CharVector{N,L}(Vector{UInt8}(undef, bufferlength), 1, N)
end

#function Base.similar(cv::CharVector{N,L}, nrow::Integer) where {N,L}
#    bufferlength = nrow * N
#    return CharVector{N,nrow}(Vector{UInt8}(undef, bufferlength), 1, N)
#end

function Base.copy(cv::CharVector{N,L}) where {N,L}
    bufferlength = L * N
    buffer = Vector{UInt8}(undef, bufferlength)
    offset = cv.offset
    recordlength = cv.recordlength
    for i = 0:(L-1)
        buffer[i*N+1:i*N+N] = cv.buffer[i*recordlength+offset:i * recordlength+offset+N-1]
    end
    return CharVector{N,L}(buffer, 1, N)
end

Tables.allocatecolumn(::Type{FixedSizeString{N}}, L) where {N} =
    CharVector{N,L}(Vector{UInt8}(undef, L * N), 1, N)
function Base.copyto!(
    dest::CharVector{N,L},
    d_o::Integer,
    src::CharVector{N,M},
) where {N,L,M}
    for i = 0:(M-1)
        dest.buffer[(d_o-1+i)*dest.recordlength+dest.offset:(d_o - 1 + i) *
                                                            dest.recordlength+dest.offset+N-1] =
            src.buffer[i*src.recordlength+src.offset:i * src.recordlength+src.offset+N-1]
    end
end

# TODO: deleteat!(cv::CharVector, i::Integer)
# TODO: deleteat!(cv::CharVector, inds)

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
types(t::FWFTable) = [var.datatype for var in specs(t)]
Tables.schema(t::FWFTable) = Tables.schema(names(t), types(t))
Base.length(t::FWFTable) = length(cols(t)[Symbol(specs(t)[1].name)])

function createcolumn(::Type{FixedSizeString}, buffer, nrow, startpos, length, recordlength)
    CharVector{length, nrow}(buffer, startpos, recordlength)
end

function createcolumn(::Type{Union{Missing, Int64}}, buffer, nrow, startpos, length, recordlength)
    column = Vector{Union{Missing, Int64}}(missing, nrow)
    select_end = startpos + length - 1
    for i = 1:nrow
        setindex!(
            column,
            bytestoint(Int64, buffer[(i-1)*recordlength+startpos:(i-1)*recordlength+select_end]),
            i,
        )
    end
    return column
end

function createcolumn(::Type{Float64}, buffer, nrow, startpos, length, recordlength)
    column = Vector{Float64}(undef, nrow)
    select_end = startpos + length - 1
    for i = 1:nrow
        setindex!(
            column,
            bytestofloat(buffer[(i-1)*recordlength+startpos:(i-1)*recordlength+select_end]),
            i,
           )
    end
    return column
end

function File(filename::String, specs::Vector{Varspec})
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
    # TODO geen regeleinde bij laatste regel opvangen
    nrow = length(buffer) ÷ recordlength
    columns = Dict{Symbol,AbstractVector}()
    for spec in specselectie
        column = createcolumn(spec.datatype, buffer, nrow, spec.startpos, spec.length, recordlength)
        columns[Symbol(spec.name)] = column
    end
    return FWFTable(specselectie, columns)
end

function bytestoint(::Type{T}, b) where {T<:Integer}
    if length(b) == 0
        return missing
    end
    res::T = 0  # "     " -> 0?
    sign = false
    for c in b
        if c == 0x2d && !sign 
            sign = true
            continue
        elseif c == 0x20
            continue
        elseif (c < 0x30) || (c > 0x39)
            return missing
        end
        res = res * 10 + c - 48
    end
    if sign
      return -res
    end
    return res
end

function bytestofloat(b, dec = 0x2e)
    if length(b) == 0
        return NaN64
    end
    res::Float64 = 0
    sign = false
    pre = true
    factor = 0.1
    for c in b
    	if c == 0x2d && !sign
            sign = true
            continue
        elseif c == 0x20
            continue
        elseif c == dec
            pre = false
            continue
        elseif (c < 0x30) || (c > 0x39)
            return NaN64
        elseif pre
            res = res * 10 + c - 48
        else
            res = res + factor * (c - 48)
            factor = factor / 10
        end
    end
    if sign
      return -res
    end
    return res
end

function maakvec(buffer, nrow::Integer, start::Integer, lengte::Integer, recordlengte::Integer)
	buf = IOBuffer(buffer)
	res = Vector{FixedSizeString{lengte}}(undef, nrow)
	pos=start
	for i in 1:nrow
	    unsafe_copyto!(pointer(res, i), convert(Ptr{FixedSizeString{10}}, pointer(buffer, pos)), 1)
		pos = pos + recordlengte
	end
	return res
end

end # module

# vim: ts=4:sw=4:textwidth=89:fileencoding=utf-8
