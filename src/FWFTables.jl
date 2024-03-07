module FWFTables

import Tables
import Parsers
import Format
import FixedSizeStrings.FixedSizeString

export readbla, Varspec, FWFTable, File, write

"""
    Varspec

Struct to hold all relevant information for one variable. 
Fields are

- name: name of the variable
- datatype: Julia type of variable
- startpos: starting position of variabele in record
- length: length in bytes of variable
- decimals: for float-type the number of decimals

Predefined usable types are `FixedSizeString{N}` for string variables,
`Union{Missing,Int64}` and `IntXX` for integer variables, `FloatXX` for real numbers and
`Nothing` for variables that must be skipped. If you want to use another type for input,
then you must define `fw_convert` to convert the byte-string to the given type.

To write to a fixed width file, generator functions for the used types must be available.
Predefined function are defined for `FixedSizeString{N}`, `Union{Missing,Int64}`,
`Float64` and `Nothing`. For other functions the function `makewrite` must be defined.
"""
struct Varspec
    name::String
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
            datatype = convdatatype[lowercase(var[:type])]
            if var[:array] == ""
                push!(spec, Varspec(name, datatype, startpos, len, decimals))
                startpos += len
            else
                # vector of vars
                for varnr = parse(Int64, var[:astart]):parse(Int64, var[a:eind])
                    tmpname = name * "_" * string(varnr)
                    push!(spec, Varspec(tmpname, datatype, startpos, len, decimals))
                    startpos += len
                end
            end
        end
    end
    return spec
end


struct CharVector{N,L} <: AbstractVector{FixedSizeString{N}}
    buffer::Vector{UInt8}
    offset::Int64
    recordlength::Int64
end

Base.IndexStyle(cv::CharVector{N,L}) where {N,L} = IndexLinear()
Base.size(cv::CharVector{N,L}) where {N,L} = (L, )

function Base.getindex(cv::CharVector{N,L}, i::Integer) where {N,L}
    @boundscheck 1 <= i <= L || throw("index $i out of bounds")
    startpos = (i - 1) * cv.recordlength + cv.offset
    endpos = (i - 1) * cv.recordlength + cv.offset + N - 1
    @inbounds s = FixedSizeString{N}(view(cv.buffer, startpos:endpos))
    return s
end

#Base.view(cv::CharVector{N, L}, i::Integer) where {N, L} = cv[i]
#Base.view(cv::CharVector{N, L}, inds::AbstractUnitRange) where {N, L} = view(cv, collect(inds))
#Base.getindex(cv::CharVector{N, L}, inds::AbstractUnitRange) where {N, L} = view(cv, collect(inds))
#Base.getindex(cv::CharVector{N, L}, inds) where {N, L} = view(cv, inds)
#Base.Broadcast.dotview(cv::CharVector{N, L}, inds::AbstractUnitRange) where {N, L} = Base.Broadcast.dotview(cv, collect(inds))

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

fw_convert(T::Type{<:Integer}, value) = Parsers.parse(T, value)
fw_convert(::Type{Union{Missing, Int64}}, value) = something(Parsers.tryparse(Int64, value), missing)
fw_convert(::Type{Float64}, value) = something(Parsers.tryparse(Float64, value), NaN64)
fw_convert(::Type{Float32}, value) = something(Parsers.tryparse(Float32, value), NaN32)
fw_convert(::Type{Float16}, value) = something(Parsers.tryparse(Float16, value), NaN32)

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
        File(io, specs)
    end
end

"""
    File(io:IO, specs::Vector{Varspec})

Read a fixed width file.
"""
function File(io::IO, specs::Vector{Varspec})
    recordlength = maximum(spec.startpos + spec.length - 1 for spec in specs)
    specselectie = [x for x in specs if x.datatype !== Nothing]
    pos = position(io)
    bom = read(io, 3) == [0xef, 0xbb, 0xbf]
    if !bom
        seek(io, pos)
    end
    buffer = read(io)
    rawlength = recordlength
    while (length(buffer) > recordlength) && (buffer[recordlength+1] in [0x0a, 0x0d])
        recordlength = recordlength + 1
    end
    crlflength = recordlength - rawlength
    # geen regeleinde bij laatste regel opvangen
    if crlflength > 0 && !(buffer[end] in [0x0a, 0x0d])
        nrow = (length(buffer) + crlflength) ÷ recordlength
    else
        nrow = length(buffer) ÷ recordlength
    end
    columns = Dict{Symbol,AbstractVector}()
    Threads.@threads for spec in specselectie
        column = createcolumn(spec.datatype, buffer, nrow, spec.startpos, spec.length, recordlength)
        columns[Symbol(spec.name)] = column
    end
    return FWFTable(specselectie, columns)
end


makewrite(::Type{FixedSizeString}, spec) = (io, val) -> Base.write(io, Ref(val))
function makewrite(::Type{Union{Missing, Int64}}, spec)
    fe = Format.FormatSpec(string("0>", spec.length, "d"))
    spaces = repeat(" ", spec.length) 
    return (io, val) -> ismissing(val) ? spaces : Format.printfmt(io, fe, val)
end
function makewrite(::Type{Float64}, spec)
    fe = Format.FormatSpec(string("0>", spec.length, ".", spec.decimals, "f"))
    return (io, val) -> Format.printfmt(io, fe, val)
end
function makewrite(::Type{Nothing}, spec)
    spaces = repeat(" ", spec.length) 
    return io -> Base.write(io, spaces)
end


"""
    write(filename::String, blafilename::String, table)

Write the table to a fixed width ascii file. The blafile contains the specifications of
the file in the Blaise format.

Using the datatype, length and if relevant the decimals from the specification a conversion
function is generated to write a single value using the required width.
For String, Integer, Real and Dummy there are conversion-functions generated by
*makewrite*.  For other datatypes, the required generator must be created. The signature
is (T::Type, spec::Varspec), and the return value must be a function that accepts
(io::IO, value::T) and writes the ascii-representation of the value to io.
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
    writefcies = [Symbol(spec.name) => makewrite(spec.datatype, spec) for spec in specs]
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
