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

Base.tryparse(::Type{Int}, ::Nothing) = 0

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
            decimals = Base.tryparse(Int64, var[:decimals])
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

function makewrite(spec)
    if spec.datatype == FixedSizeString
        return (io, val) -> write(io, Ref(val))
    elseif spec.datatype == Union{Missing, Int64}
        fs = "0>" * string(spec.length) * "d"
        fe = Formatting.FormatSpec(fs)
        return (io, val) -> Formatting.printfmt(io, fe, val)
    elseif spec.datatype == Float64
        fs = "0>" * string(spec.length) * "." * String(spec.decimals) * "f"
        fe = Formatting.FormatSpec(fs)
        return (io, val) -> Formatting.printfmt(io, fe, val)
    else
        spaces = string([" " for i in 1:spec.length]) 
        return (io, val) -> write(io, spaces)
    end
end


function write(filename::String, specs::Vector{Varspec}, table)
    fe = Formatting.FormatExpr(makefmt(specs))
    writefcies = [spec.name => makewrite(spec) for spec in specs]
    open(filename, "w") do io
        for row in Tables.rows(table)
            for (name, writefcie) in writefcies
                #if name == :dummy
                writefcie(io, row[name])
            end    
            println(io, "")
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
    return Vector{FixedSizeString{N}}(undef, L)
end

function Base.copy(cv::CharVector{N,L}) where {N,L}
	res = Vector{FixedSizeString{N}}(undef, L)
    recordlength = cv.recordlength
    pointer_from = pointer(cv.buffer, 1) + cv.offset - 1
    pointer_to = convert(Ptr{UInt8}, pointer(res, 1))
	for i in 1:L
        unsafe_copyto!(pointer_to, pointer_from, N)
        pointer_from = pointer_from + recordlength
        pointer_to = pointer_to + N
	end
	return res
end

Tables.allocatecolumn(::Type{FixedSizeString{N}}, L) where {N} =
    Vector{FixedSizeString{N}}(undef, L)

function Base.copyto!(
    dest::Vector{FixedSizeString{N}},
    d_o::Integer,
    src::CharVector{N,M},
) where {N,M}
    recordlength = src.recordlength
    pointer_from = pointer(src.buffer, 1) + src.start - 1
    pointer_to = convert(Ptr{UInt8}, pointer(dest, d_o))
    for i = 1:M
        unsafe_copyto!(pointer_to, pointer_from, N)
        pointer_from = pointer_from + recordlength
        pointer_to = pointer_to + N
    end
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
types(t::FWFTable) = [var.datatype for var in specs(t)]
Tables.schema(t::FWFTable) = Tables.schema(names(t), types(t))
Base.length(t::FWFTable) = length(cols(t)[Symbol(specs(t)[1].name)])

function createcolumn(::Type{FixedSizeString}, buffer, nrow, startpos, length, recordlength)
    CharVector{length, nrow}(buffer, startpos, recordlength)
end

fw_convert(::Type{Union{Missing, Int64}}, value) = something(Parsers.tryparse(Int64, value), missing)
fw_convert(::Type{Float64}, value) = something(Parsers.tryparse(Float64, value), NaN64)

function createcolumn(T, buffer, nrow, startpos, length, recordlength)
    column = Vector{T}(undef, nrow)
    pos = startpos
    endpos = startpos + length - 1
    for i = 1:nrow
        setindex!(column, fw_convert(T, view(buffer, pos:endpos)), i)
        pos = pos + recordlength
        endpos = endpos + recordlength
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


function maakvec2(buffer, nrow::Integer, start::Integer, lengte::Integer, recordlengte::Integer)
	res = Vector{FixedSizeString{lengte}}(undef, nrow)
    pointer_from = pointer(buffer, 1) + start - 1
    pointer_to = convert(Ptr{UInt8}, pointer(res, 1))
	for i in 1:nrow
        unsafe_copyto!(pointer_to, pointer_from, lengte)
        pointer_from = pointer_from + recordlengte
        pointer_to = pointer_to + lengte
	end
	return res
end

function maakvec(buffer, nrow::Integer, start::Integer, lengte::Integer, recordlengte::Integer)
	res = Vector{FixedSizeString{lengte}}(undef, nrow)
    pointer_from = convert(Ptr{FixedSizeString{lengte}}, pointer(buffer, 1)) + start - 1
    pointer_to = pointer(res, 1)
	for i in 1:nrow
        unsafe_copyto!(pointer_to, pointer_from, 1)
        pointer_from = pointer_from + recordlengte
        pointer_to = pointer_to + lengte
	end
	return res
end

end # module

# vim: ts=4:sw=4:textwidth=89:fileencoding=utf-8
