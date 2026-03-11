module AstroUnitFormats

export parse_unit, unit_string, CDS, VOUnit, FITS

using Unitful: Unitful, NoUnits, uparse
using UnitfulAstro: UnitfulAstro
using UnitfulAngles: UnitfulAngles

const _UNIT_CONTEXT = [UnitfulAngles, UnitfulAstro, Unitful]


# ── Name mapping dicts ────────────────────────────────────────────────────────
# Each format maps format-specific unit names → Unitful-parseable names.

# Shared base: names that need remapping in all formats
const _BASE_NAME_DICT = Dict(
    "Angstroem" => "Å",    "Angstrom" => "Å",    "angstrom" => "Å",    "AA" => "Å",
    "solMass" => "Msun",    "solRad" => "Rsun",    "solLum" => "Lsun",
    "arcsec" => "arcsecond","arcmin" => "arcminute",
    "Ohm" => "Ω",
    "lyr" => "ly",
    "au" => "AU",
    "um" => "μm",
    "min" => "minute",      "h" => "hr",
)

# CDS-only additions
const _CDS_NAME_DICT = merge(_BASE_NAME_DICT, Dict(
    "jovMass" => "Mjup",    "geoMass" => "Mearth", "Mgeo" => "Mearth", "Rgeo" => "Rearth",
    "arcs" => "arcsecond",  "arcm" => "arcminute",
    "uarcsec" => "μas",     "marcsec" => "mas",
    "gauss" => "Gauss",
    "al" => "ly",
    "eps0" => "ε0",         "mu0" => "μ0",
    "a0" => "bohrRadius",
    "sec" => "s",
))

# Reverse mapping: Unitful abbreviation → canonical CDS name (only where they differ)
const _CDS_REVERSE_NAME_DICT = Dict(
    "M⊙" => "solMass",  "R⊙" => "solRad",  "L⊙" => "solLum",
    "M⊕" => "Mearth",   "R⊕" => "Rearth",
    "Å" => "Angstrom",   "Ω" => "Ohm",      "°" => "deg",
    "hr" => "h",          "minute" => "min",
    "′" => "arcmin",      "″" => "arcsec",
    "ly" => "lyr",
)

# VOUnit additions (IVOA VOUnits 1.1)
const _VOUNIT_NAME_DICT = merge(_BASE_NAME_DICT, Dict(
    "G" => "Gauss",         # VOUnit G = Gauss (deprecated); CDS G = gravitational constant (no remap needed)
    "a" => "yr",            # Julian year — Unitful a = Are (area)
))

# FITS additions (FITS Standard 4.0, Section 4.3)
const _FITS_NAME_DICT = merge(_BASE_NAME_DICT, Dict(
    "G" => "Gauss",         # same as VOUnit
    "a" => "yr",
))


# ── Regex helpers ─────────────────────────────────────────────────────────────

# Letter-only boundaries (not \b) since \b treats digits as word chars,
# so e.g. "10pix" has no \b between "10" and "pix".
# Alternation sorted longest-first so e.g. "arcsec" matches before "arcs".
_letter_boundary_re(alts) = Regex("(?<![a-zA-Z])(" * join(sort!(collect(alts), by=length, rev=true), "|") * ")(?![a-zA-Z])")

const _CDS_NAME_RE    = _letter_boundary_re(keys(_CDS_NAME_DICT))
const _VOUNIT_NAME_RE = _letter_boundary_re(keys(_VOUNIT_NAME_DICT))
const _FITS_NAME_RE   = _letter_boundary_re(keys(_FITS_NAME_DICT))

# Generic: union of all format dicts; VOUnit wins on conflicts (G=Gauss, a=yr)
const _GENERIC_NAME_DICT = merge(_CDS_NAME_DICT, _VOUNIT_NAME_DICT)
const _GENERIC_NAME_RE   = _letter_boundary_re(keys(_GENERIC_NAME_DICT))


# ── Dimensionless unit sets ───────────────────────────────────────────────────
# Same across all formats: units with no Unitful equivalent → treated as "1".

const _DIMENSIONLESS = Set([
    "pix", "pixel", "ct", "count", "photon", "ph",
    "adu", "DN", "dn", "electron",
    "chan", "bin", "beam", "vox", "voxel",
    "bit", "byte",
])
# Matches a dimensionless unit + optional trailing power ([+-]?\d+),
# so that e.g. "beam-1" is replaced as a whole instead of leaving an orphaned "-1".
_sorted_dimless = sort!(collect(_DIMENSIONLESS), by=length, rev=true)
_dimless_alt = join(_sorted_dimless, "|")
# After a digit: "10pix" → keep digit, drop unit (10×1=10). Avoids "10"*"1"="101".
# (?!\.) in the power group prevents consuming digits before a decimal point (e.g. "pix0.1nm").
# (?![\da-zA-Z]) prevents matching when the unit is directly followed by a digit.
const _DIMLESS_AFTER_DIGIT = Regex("(\\d)(" * _dimless_alt * ")([+-]?\\d+(?!\\.))?(?![\\da-zA-Z])")
# Otherwise: replace with "1".
const _DIMLESS_RE = Regex("(?<![a-zA-Z])(" * _dimless_alt * ")([+-]?\\d+(?!\\.))?(?![\\da-zA-Z])")


# ── Shared pipeline ──────────────────────────────────────────────────────────
# All syntax transformations stacked; format-irrelevant ones are no-ops.

function _parse_unit(s::AbstractString, name_re, name_dict)
    s = strip(s)

    # VOUnit quoted unit strings: strip single quotes (e.g. '1' → 1, 'dex' → dex)
    if occursin('\'', s)
        @warn "Stripping single quotes from unit string (VOUnit quoted-unit notation)" unit_string=s
        s = replace(s, "'" => "")
    end

    # Log notation: [unit] (CDS brackets) or log(unit) → exp10 transform
    m = match(r"^(?:\[(.+)\]|log\((.+)\))$", s)
    if m !== nothing
        s = string(something(m[1], m[2]))
        valuefn = exp10
    else
        valuefn = identity
    end

    # Dimensionless markers (CDS: "---"/"-", VOUnit/FITS: "unknown"/"UNKNOWN", empty)
    if s == "---" || s == "-" || s == "unknown" || s == "UNKNOWN" || isempty(s)
        return (; unit = NoUnits, valuefn)
    end

    # Percent: bare "%" or "20%" etc.
    s = replace(s, "%" => "*percent")
    startswith(s, '*') && (s = s[2:end])

    # CDS backslash-escaped constants: \h = Planck constant (must be before name dict)
    # PLANCKh placeholder has letters on both sides of "h", so the letter-boundary
    # regex won't match it during the name-dict pass.
    s = replace(s, r"\\h(?![a-zA-Z])" => "PLANCKh")

    # Format-specific name → Unitful name (single-pass via precompiled alternation regex)
    s = replace(s, name_re => m -> name_dict[m])

    # Restore Planck constant: PLANCKh → h (Unitful symbol for Planck constant)
    s = replace(s, "PLANCKh" => "h")

    # Factor × 10^exp: "1.5x10+11m" → "1.5e+11*m", "1.5×10+11/m" → "1.5e+11/m"
    m = match(r"^([\d.]+)[x×]10([+-]\d+)(.*)", s)
    if m !== nothing
        factor, exp, rest = string.(m.captures)
        sep = isempty(rest) || startswith(rest, '/') ? "" : "*"
        s = "$(factor)e$(exp)$(sep)$rest"
    else
        # Pure 10^exp prefix: "10+20cm-2" → "1e+20*cm-2", "10+20/m" → "1e+20/m"
        m = match(r"^10([+-]\d+)(.*)", s)
        if m !== nothing
            exp, rest = string.(m.captures)
            sep = isempty(rest) || startswith(rest, '/') ? "" : "*"
            s = "1e$(exp)$(sep)$rest"
        end
    end

    # Leading division: "/s" → "s^-1"
    startswith(s, "/") && (s = s[2:end] * "^-1")

    # ** → ^ (VOUnit/FITS power syntax; no-op for CDS which has no **)
    s = replace(s, "**" => "^")

    # Dot multiplication (only when letter on at least one side, preserving decimals)
    s = replace(s, r"\.(?=[a-zA-Z])" => "*")
    s = replace(s, r"(?<=[a-zA-Z])\." => "*")

    # Dimensionless units → replaced by 1 (pix, ct, photon, adu, etc.)
    # Done after dot→* so that e.g. "Jy.beam-1" is already "Jy*beam-1" and won't produce "1.1".
    s = replace(s, _DIMLESS_AFTER_DIGIT => s"\1")  # "10pix" → "10" (10×1=10)
    s = replace(s, _DIMLESS_RE => "1")              # "beam-1" → "1"

    # Bare-sign powers: "s-2" → "s^-2" (CDS; after ** already converted, this is a no-op for VOUnit/FITS)
    # Negative lookbehind for digit to skip 'e' in scientific notation like "1e+21"
    s = replace(s, r"(?<!\d)([a-zA-Z])([+-]\d+)" => s"\1^\2")

    # Unsigned powers: "m2" → "m^2" (letter+digit at end or before operator)
    s = replace(s, r"([a-zA-Z])(\d+)(?=$|[/*.)])" => s"\1^\2")

    try
        u = uparse(s; unit_context=_UNIT_CONTEXT)
        return (; unit = u, valuefn)
    catch
        @warn "Could not parse unit string" unit_string=s
        return (; unit = NoUnits, valuefn)
    end
end


# ── Format types ──────────────────────────────────────────────────────────────

struct CDS end
struct VOUnit end
struct FITS end


# ── Public API ────────────────────────────────────────────────────────────────

"""
    parse_unit(s) -> (; unit, valuefn)
    parse_unit(s, CDS()) -> (; unit, valuefn)
    parse_unit(s, VOUnit()) -> (; unit, valuefn)
    parse_unit(s, FITS()) -> (; unit, valuefn)

Parse an astronomical unit string. Without a format argument, uses a generic
parser that accepts all syntax variants (CDS, VOUnit, FITS).

Returns a `NamedTuple` with:
- `unit`: a `Unitful.FreeUnits` object (`Unitful.NoUnits` for dimensionless or unparseable)
- `valuefn`: `identity` for regular units, `exp10` for log/bracket units like `[K]` or `log(Hz)`

# Examples
```julia
parse_unit("km/s")                       # (unit = km s⁻¹, valuefn = identity)
parse_unit("km.s-1", CDS())              # CDS dot-multiplication, bare-sign powers
parse_unit("m**-2", VOUnit())            # VOUnit ** power syntax
parse_unit("erg.s**-1.cm**-2", FITS())   # FITS format
parse_unit("[K]")                         # (unit = K, valuefn = exp10)
```
"""
parse_unit(s::AbstractString) = _parse_unit(s, _GENERIC_NAME_RE, _GENERIC_NAME_DICT)
parse_unit(s::AbstractString, ::CDS) = _parse_unit(s, _CDS_NAME_RE, _CDS_NAME_DICT)
parse_unit(s::AbstractString, ::VOUnit) = _parse_unit(s, _VOUNIT_NAME_RE, _VOUNIT_NAME_DICT)
parse_unit(s::AbstractString, ::FITS) = _parse_unit(s, _FITS_NAME_RE, _FITS_NAME_DICT)


"""
    unit_string(unit) -> String

Convert a `Unitful.FreeUnits` to a CDS-format unit string.
Returns `"---"` for dimensionless (`NoUnits`).

# Examples
```julia
unit_string(u"km/s")    # "km.s**-1"
unit_string(u"Msun")   # "solMass"
unit_string(NoUnits)   # "---"
```
"""
function unit_string(unit::Unitful.Units)
    comps = typeof(unit).parameters[1]
    isempty(comps) && return "---"
    join(map(comps) do comp
        name = get(_CDS_REVERSE_NAME_DICT, Unitful.abbr(comp), Unitful.abbr(comp))
        s = Unitful.prefix(comp) * name
        p = Int(Unitful.power(comp))
        p == 1 ? s : s * "**" * string(p)
    end, ".")
end

function unit_string(unit::Unitful.MixedUnits)
    unit.units == NoUnits || error("unit_string does not support compound MixedUnits like $unit")
    string(unit)
end

end
