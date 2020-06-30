module FWFTables

import Base.iterate
import Tables
import Parsers
import Formatting

export readlba, Blavar, FWFTable, makefmt, File, write


struct Blavar
  name::String
  slice::UnitRange{Int64}
  datatype::DataType
  length::Int64
  decimals::Int64
end

varregex = r"(?i)\s*((?P<name>\w+)\s*(?P<tekst>\".*\"|)\s*:|)\s*(?P<array>ARRAY\[(?P<astart>\d+)\.\.(?P<aeind>\d+)\]\s*OF\s*|)(?P<type>DUMMY|STRING|REAL|INTEGER)\s*\[\s*(?P<length>\d+)(\s*,\s*(?P<decimals>\d+)|)\s*\].*"

convdatatype = Dict(
  "string" => String,
  "integer" => Int64,
  "dummy" => Nothing,
  "real" => Float64,
)


Parsers.tryparse(::Type{Int}, ::Nothing) = 0
Parsers.tryparse(::Type{String}, s::String) = s


function readbla(filename)
  bla = Vector{Blavar}()
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
      datatype = convdatatype[lowercase(var[:type])]
      len = parse(Int64, var[:length])
      decimals = Parsers.tryparse(Int64, var[:decimals])
      if var[:array] == ""
        slice = startpos:(startpos+len-1)
        push!(bla, Blavar(name, slice, datatype, len, decimals))
        startpos += len
      else
        # vector of vars
        for varnr = parse(Int64, var[:astart]):parse(Int64, var[a:eind])
          tmpname = name * "_" * string(varnr)
          slice = startpos:(startpos+len-1)
          push!(bla, Blavar(tmpname, slice, datatype, len, decimals))
          startpos += len
        end
      end
    end
  end
  return bla
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
end

function makefmt(bla::Vector{Blavar})
  string([makefmt(x) for x in bla]...)
end

struct FWFTable
  handle::IOStream
  bla::Vector{Blavar}
  names::Dict{Symbol,Int}
  numberofrecords::Int
  recordlength::Int
  crlflength::Int
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
names(f::FWFTable) = [var.name for var in f.bla]
Tables.columnnames(f::FWFTable) = [var.name for var in f.bla]
types(f::FWFTable) = [var.datatype for var in f.bla]
Tables.schema(f::FWFTable) = Tables.Schema(names(f), types(f))


Tables.rowaccess(::FWFTable) = true
Tables.rows(f::FWFTable) = f
Base.eltype(f::FWFTable) = FWFTableRow
Base.length(f::FWFTable) = getfield(f, :numberofrecords)
Base.iterate(f::FWFTable, st = 0) = nextline(f, st)

Tables.getcolumn(r::FWFTableRow, s::String) = Tables.getcolumn(r, Symbol(s))

Tables.columnnames(r::FWFTableRow) = names(getfield(r, :source))

function Tables.getcolumn(r::FWFTableRow, col::Int)
  var = getfield(getfield(r, :source), :bla)[col]
  value = Parsers.tryparse(var.datatype, getfield(r, :rawrow)[var.slice])
  if isnothing(value)
    return missing
  else
    return value
  end
end

function Tables.getcolumn(r::FWFTableRow, nm::Symbol)
  col = getfield(getfield(r, :source), :names)[nm]
  var = getfield(getfield(r, :source), :bla)[col]
  rawrow = getfield(r, :rawrow)
  if var.slice.stop <= length(rawrow)
    value = Parsers.tryparse(var.datatype, rawrow[var.slice])
  else
    value = Nothing
  end
  if isnothing(value)
    return missing
  else
    return value
  end
end

function File(filename::String, blafilename::String, crlf = 1)
  bla = readbla(blafilename)
  File(filename, bla, crlf)
end

function File(filename::String, bla::Vector{Blavar}, crlf = 1)
  blaselectie = [x for x in bla if x.datatype !== Nothing]
  d = Dict([Symbol(elt.name) => idx for (idx, elt) in enumerate(blaselectie)])
  recordlength = sum(elt.length for elt in bla)
  filesize = stat(filename).size
  if filesize % (recordlength + crlf) in (0, crlf)
    numberofrecords = filesize ÷ (recordlength + crlf)
  elseif (filesize + crlf) % (recordlength + crlf) == 0
    numberofrecords = (filesize + crlf) ÷ (recordlength + crlf)
  else
    numberofrecords = -1
  end
  handle = open(filename)
  return FWFTable(handle, blaselectie, d, numberofrecords, recordlength, crlf)
end

function write(filename::String, blafilename::String, table)
  bla::Vector{Blavar} = readbla(blafilename)
  write(filename, bla, table)
end

function write(filename::String, bla::Vector{Blavar}, table)
  fe = Formatting.FormatExpr(makefmt(bla))
  open(filename, "w") do io
    for row in Tables.rows(table)
      Formatting.printfmtln(io, fe, row...)
    end
  end
end

end # module

