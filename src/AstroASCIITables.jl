module AstroASCIITables

using StructArrays
import FixedWidthTables as FWT


# ── Column spec ──────────────────────────────────────────────────────────────

struct ColSpec
    name::Symbol
    range::UnitRange{Int}             # 1-based byte range, as in CDS format
    base_type::Type                   # Int64, Float64, or String
    nullval::Union{String,Nothing}    # nothing = not nullable; "" = blank; "-" = dashes; else exact
    description::String
    units::String
end


# ── Section delimiter (rows of only "-" or "=") ──────────────────────────────

_is_section_delimiter(line) =
    (startswith(line, "------") || startswith(line, "=======")) && allequal(strip(line))


# ── CDS format code → Julia type ─────────────────────────────────────────────

function _cds_type(fmt)
    c = lowercase(first(match(r"[a-zA-Z]", fmt).match))
    c == 'i' && return Int64
    c in ('f', 'e') && return Float64
    return String
end


# ── Simple fnmatch (only * wildcard needed for CDS filenames) ────────────────

_cds_fnmatch(pattern, name) =
    occursin(Regex("^" * replace(pattern, "." => "\\.", "*" => ".*", "?" => ".") * "\$"), name)


# ── Parse nullable marker from a column description ──────────────────────────
# CDS uses "?" to mark optional columns; optionally followed by "=<nullval>"
# and an order specifier (+, -, +=, -=).
# Returns (nullval, cleaned_description).

const _NULLABLE_RE = r"^(?:[\[\]]\S*[\[\]])?\?(?:=(\S*))?[-+]?=?\s*(\S.*)?$"

function _parse_nullable(descr)
    m = match(_NULLABLE_RE, descr)
    m === nothing && return (nothing, descr)  # not nullable
    nullval = m[1] === nothing ? "" :         # plain "?": blank field means missing
              m[1] == "-"     ? "-" :         # "?=-": any run of dashes means missing
                                string(m[1])  # "?=-9.9" etc: exact sentinel value
    (nullval, m[2] === nothing ? "" : string(strip(m[2])))
end


# ── Column-definition line regex (ported from astropy cds.py re_col_def) ─────

const _COL_RE = r"^\s*(?:(?<start>\d+)\s*-)?\s*(?<stop>\d+)\s+(?<fmt>[\w.]+)\s+(?<units>\S+)\s+(?<label>\S+)(?:\s+(?<descr>\S.*?))?$"

function _try_parse_col_line(line)
    m = match(_COL_RE, line)
    m === nothing && return nothing
    stop   = parse(Int, m[:stop])
    start  = m[:start] === nothing ? stop : parse(Int, match(r"\d+", m[:start]).match)
    nullval, descr = _parse_nullable(string(something(m[:descr], "")))
    ColSpec(Symbol(m[:label]), start:stop, _cds_type(m[:fmt]),
            nullval, string(strip(descr)), string(m[:units]))
end


# ── Extract the header section for a specific table from a ReadMe ────────────
# Mirrors astropy's get_cols ReadMe-mode: scan for matching Byte-by-byte header,
# collect lines until the 3rd section delimiter.

function _extract_readme_section(lines, tablename)
    result = String[]
    in_header = false
    delim_count = 0
    for line in lines
        if in_header
            push!(result, line)
            if _is_section_delimiter(line)
                delim_count += 1
                delim_count == 3 && break
            end
        else
            m = match(r"Byte-by-byte Description of file:\s*(.+)"i, line)
            if m !== nothing
                # File list may be comma/space separated and may use "*" globs
                names = filter(!isempty, split(strip(m[1]), r"[,\s]+"))
                if any(pat -> _cds_fnmatch(string(pat), tablename), names)
                    in_header = true
                    push!(result, line)
                end
            end
        end
    end
    isempty(result) && error("Can't find table \"$tablename\" in ReadMe")
    return result
end


# ── Parse the Byte-by-byte Description block into ColSpecs ───────────────────

_append_descr(col::ColSpec, text) =
    ColSpec(col.name, col.range, col.base_type, col.nullval,
            col.description * (isempty(col.description) ? "" : " ") * text,
            col.units)

function _parse_header(lines, tablename=nothing)
    work = tablename === nothing ? lines : _extract_readme_section(lines, tablename)

    # Find the index of the last consecutive "Byte-by-byte Description" line.
    # Some ReadMe files stack two Byte-by-byte lines for the same column block.
    i_bbb = nothing
    for (i, line) in enumerate(work)
        if occursin(r"Byte-by-byte Description"i, line)
            i_bbb = i
        elseif i_bbb !== nothing
            break  # first non-Byte-by-byte line after the block
        end
    end
    i_bbb === nothing && error("No \"Byte-by-byte Description\" found")

    # Layout: i_bbb → last Byte-by-byte line
    #   +1: section delimiter (-------)
    #   +2: "Bytes Format Units Label Explanations" header row
    #   +3: section delimiter (-------)
    #   +4: first column definition line
    col_lines = Iterators.takewhile(!_is_section_delimiter, work[i_bbb+4:end])
    cols = foldl(col_lines; init=ColSpec[]) do acc, line
        col = _try_parse_col_line(line)
        if col !== nothing
            push!(acc, col)
        elseif !isempty(acc) && !isempty(strip(line))
            # Continuation line: append to previous column's description
            acc[end] = _append_descr(acc[end], strip(line))
        end
        acc
    end

    isempty(cols) && error("No columns found in header")
    return cols
end


# ── Build FixedWidthTables colspecs from ColSpec list ────────────────────────

_fwt_specs(cols) = NamedTuple{Tuple(c.name for c in cols)}(
    Tuple((c.range, !isnothing(c.nullval) ? String : c.base_type) for c in cols))


# ── Stream header lines from a self-contained file, leave IO at data start ───
# Reads lines one at a time, stops (and seeks back) when the first data row is
# reached.  Data starts after the 3rd section delimiter inside the BBB block.

function _collect_header_streaming!(io::IO)
    header = String[]
    in_bbb = false
    bbb_delims = 0
    in_notes = false
    last_delim_end_pos = nothing   # byte position right after the last post-coldef delimiter
    while !eof(io)
        pos = position(io)
        line = readline(io)
        if !in_bbb || bbb_delims < 3
            # Preamble or inside the Byte-by-byte Description block
            push!(header, line)
            if occursin(r"Byte-by-byte Description"i, line)
                in_bbb, bbb_delims = true, 0
            elseif in_bbb && _is_section_delimiter(line)
                bbb_delims += 1
            end
        else
            # Post-coldef: anything here is Notes or data
            if _is_section_delimiter(line)
                push!(header, line)
                in_notes = false
                last_delim_end_pos = position(io)   # start of next line
            elseif in_notes || startswith(line, "Note") || isempty(strip(line))
                push!(header, line)
                startswith(line, "Note") && (in_notes = true)
            else
                # First non-Note/non-blank line after the column block → data
                seek(io, something(last_delim_end_pos, pos))
                break
            end
        end
    end
    return _parse_header(header, nothing)
end


# ── Post-process a nullable column (read as String, convert to target type) ──

function _postprocess_col(vals, spec::ColSpec)
    T  = spec.base_type
    nv = spec.nullval
    is_null(s) = nv == ""  ? isempty(s) :
                 nv == "-" ? (!isempty(s) && all(==('-'), s)) :
                             s == nv
    map(vals) do v
        ismissing(v) && return missing
        s = strip(v)
        is_null(s) ? missing : T == String ? string(s) : parse(T, s)
    end
end


# ── Public API ────────────────────────────────────────────────────────────────

"""
    read_cds(path; readme=nothing) -> StructArray

Read a CDS/Vizier fixed-width ASCII table.

- `readme=nothing`: self-contained file (header + data in the same file).
- `readme=<path>`: ReadMe file supplies the column header; `path` is the data file.
"""
function read_cds(path; readme=nothing)
    raw, cols = if readme !== nothing
        # Header from ReadMe (lazy iteration); data file streamed separately.
        cols = open(string(readme)) do io
            _parse_header(eachline(io), basename(string(path)))
        end
        open(string(path), "r") do io
            (FWT.read(io, StructArray, _fwt_specs(cols);
                      missings=[""], allow_shorter_lines=true), cols)
        end
    else
        # Self-contained: stream header, then hand the same IO to FWT for data.
        open(string(path), "r") do io
            cols = _collect_header_streaming!(io)   # advances io past header
            (FWT.read(io, StructArray, _fwt_specs(cols);
                      missings=[""], allow_shorter_lines=true), cols)
        end
    end

    fwt_names = Tuple(c.name for c in cols)
    result_cols = map(cols) do col
        vals = getproperty(raw, col.name)
        !isnothing(col.nullval) ? _postprocess_col(vals, col) : vals
    end
    StructArray(NamedTuple{fwt_names}(Tuple(result_cols)))
end

end
