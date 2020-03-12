module FWFTables

import Base.tryparse, Base.iterate
import Tables
import Tables.columnnames, Tables.getcolumn


struct Blavar
  name::String
  slice::UnitRange{Int64}
  datatype::DataType
  length::Int64
  decimals::Int64
end

varregex = r"(?i)\s*((?P<name>\w+)\s*(?P<tekst>\".*\"|)\s*:|)\s*(?P<array>ARRAY\[(?P<astart>\d+)\.\.(?P<aeind>\d+)\]\s*OF\s*|)(?P<type>DUMMY|STRING|REAL|INTEGER)\s*\[\s*(?P<length>\d+)(\s*,\s*(?P<decimals>\d+)|)\s*\].*"

convdatatype = Dict(
                 "string"=>String, 
                 "integer"=>Int64,
                 "dummy"=>String,
                 "real"=>Float64
                 )


tryparse(::Type{Int}, ::Nothing) = 0
tryparse(::Type{String}, s::String) = s 


function readbla(filename)
  bla = Vector{Blavar}()
  file = open(filename)
  regels = readlines(file)
  close(file)
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
      name = lowercase(var[:name])
      datatype = convdatatype[lowercase(var[:type])]
      len = parse(Int64, var[:length])
      decimals = tryparse(Int64, var[:decimals])
      slice=startpos:(startpos+len-1)
      #TODO: array naar aparte variabelen mappen
      push!(bla, Blavar(name, slice, datatype, len, decimals))
      startpos += len
    end
  end
  return bla
end 



struct FWFTable 
  handle::IOStream
  bla::Vector{Blavar}
  names::Dict{Symbol, Int}
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
Base.iterate(f::FWFTable, st=0) = nextline(f, st)

Tables.getcolumn(r::FWFTableRow, s::String) = Tables.getcolumn(r, Symbol(s))

Tables.columnnames(r::FWFTableRow) = names(getfield(r, :source))

function Tables.getcolumn(r::FWFTableRow, col::Int)
  var = getfield(getfield(r, :source), :bla)[col]
  value = tryparse(var.datatype, getfield(r, :rawrow)[var.slice])
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
    value = tryparse(var.datatype, rawrow[var.slice])
  else
    value = Nothing
  end
  if isnothing(value)
    return missing
  else
    return value
  end
end

function File(filename, blafilename, crlf=1)
  bla = readbla(blafilename)
  d = Dict([Symbol(elt.name)=>idx for (idx, elt) in enumerate(bla)])
  recordlength = sum(elt.length for elt in bla)
  filesize = stat(filename).size
  if  filesize % (recordlength + crlf) in (0, crlf)
    numberofrecords = filesize ÷ (recordlength + crlf)
  elseif (filesize + crlf) % (recordlength + crlf) == 0
    numberofrecords = (filesize + crlf) ÷ (recordlength + crlf)
  else
    numberofrecords = -1
  end
  handle = open(filename)
  return FWFTable(handle, bla, d, numberofrecords, recordlength, crlf)
end
end # module
