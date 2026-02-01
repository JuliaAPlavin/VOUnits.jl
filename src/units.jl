using Unitful: Unitful, NoUnits, uparse
using UnitfulAstro: UnitfulAstro
using UnitfulAngles: UnitfulAngles

const _CDS_UNIT_CONTEXT = [UnitfulAngles, UnitfulAstro, Unitful]

# CDS name → Unitful name: single precompiled regex + dict lookup
const _CDS_NAME_DICT = Dict(
    "Angstroem" => "Å",    "Angstrom" => "Å",    "AA" => "Å",
    "jovMass" => "Mjup",    "geoMass" => "Mearth", "Mgeo" => "Mearth", "Rgeo" => "Rearth",
    "solMass" => "Msun",    "solRad" => "Rsun",    "solLum" => "Lsun",
    "uarcsec" => "μas",     "marcsec" => "mas",       # UnitfulAngles has mas/μas as dedicated units
    "arcsec" => "arcsecond","arcmin" => "arcminute","arcs" => "arcsecond","arcm" => "arcminute",
    "gauss" => "Gauss",     "Ohm" => "Ω",
    "lyr" => "ly",          "al" => "ly",
    "eps0" => "ε0",         "mu0" => "μ0",           # Unitful physical constants (Unicode)
    "a0" => "bohrRadius",                              # Bohr radius: prevent a0 → a^0 misparsing
    "au" => "AU",                                      # CDS accepts both AU and au
    "sec" => "s",           "min" => "minute",     "h" => "hr",
)
# Letter-only boundaries (not \b) since \b treats digits as word chars,
# so e.g. "10pix" has no \b between "10" and "pix".
# Alternation sorted longest-first so e.g. "arcsec" matches before "arcs".
_letter_boundary_re(alts) = Regex("(?<![a-zA-Z])(" * join(sort!(collect(alts), by=length, rev=true), "|") * ")(?![a-zA-Z])")
const _CDS_NAME_RE = _letter_boundary_re(keys(_CDS_NAME_DICT))

# CDS units with no Unitful equivalent → treated as dimensionless ("1")
const _CDS_DIMENSIONLESS = Set([
    "pix", "pixel", "ct", "count", "photon", "ph",
    "adu", "DN", "dn", "electron",
    "chan", "bin", "beam", "vox", "voxel",
    "bit", "byte",
])
# Matches a dimensionless unit + optional trailing power ([+-]?\d+),
# so that e.g. "beam-1" is replaced as a whole instead of leaving an orphaned "-1".
_sorted_dimless = sort!(collect(_CDS_DIMENSIONLESS), by=length, rev=true)
_dimless_alt = join(_sorted_dimless, "|")
# After a digit: "10pix" → keep digit, drop unit (10×1=10). Avoids "10"*"1"="101".
const _CDS_DIMLESS_AFTER_DIGIT = Regex("(\\d)(" * _dimless_alt * ")([+-]?\\d+)?(?![a-zA-Z])")
# Otherwise: replace with "1".
const _CDS_DIMENSIONLESS_RE = Regex("(?<![a-zA-Z])(" * _dimless_alt * ")([+-]?\\d+)?(?![a-zA-Z])")

"""
    cds_to_unitful(s) -> (unit, valuefn)

Convert a CDS unit string to a Unitful unit and a value transformation function.

Returns a `NamedTuple` with:
- `unit`: a `Unitful.FreeUnits` object (`Unitful.NoUnits` for dimensionless or unparseable)
- `valuefn`: `identity` for regular units, `exp10` for dex/bracket units like `[K]`

# Examples
```julia
cds_to_unitful("km.s-1")   # (unit = km s⁻¹, valuefn = identity)
cds_to_unitful("[K]")       # (unit = K, valuefn = exp10)
cds_to_unitful("---")       # (unit = , valuefn = identity)
```
"""
function cds_to_unitful(s::AbstractString)
    s = strip(s)

    # Bracket/dex notation: [K] → K with exp10 transform
    m = match(r"^\[(.+)\]$", s)
    valuefn = m !== nothing ? exp10 : identity
    m !== nothing && (s = string(m.captures[1]))

    # Dimensionless
    (s == "---" || s == "-") && return (; unit = NoUnits, valuefn)

    # Percent: bare "%" or "20%" etc.
    s = replace(s, "%" => "*percent")
    startswith(s, '*') && (s = s[2:end])

    # CDS backslash-escaped constants: \h = Planck constant (must be before name dict)
    # PLANCKh placeholder has letters on both sides of "h", so the letter-boundary
    # regex in _CDS_NAME_RE won't match it during the name-dict pass.
    s = replace(s, r"\\h(?![a-zA-Z])" => "PLANCKh")

    # CDS name → Unitful name (single-pass via precompiled alternation regex)
    s = replace(s, _CDS_NAME_RE => m -> _CDS_NAME_DICT[m])

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

    # Dot multiplication (only when letter on at least one side, preserving decimals)
    s = replace(s, r"\.(?=[a-zA-Z])" => "*")
    s = replace(s, r"(?<=[a-zA-Z])\." => "*")

    # Known dimensionless CDS units → replaced by 1 (pix, ct, photon, adu, etc.)
    # Done after dot→* so that e.g. "Jy.beam-1" is already "Jy*beam-1" and won't produce "1.1".
    s = replace(s, _CDS_DIMLESS_AFTER_DIGIT => s"\1")  # "10pix" → "10" (10×1=10)
    s = replace(s, _CDS_DIMENSIONLESS_RE => "1")        # "beam-1" → "1"

    # Powers with sign: "s-2" → "s^-2"
    # Negative lookbehind for digit to skip 'e' in scientific notation like "1e+21"
    s = replace(s, r"(?<!\d)([a-zA-Z])([+-]\d+)" => s"\1^\2")

    # Powers without sign: "m2" → "m^2" (letter+digit at end or before operator)
    s = replace(s, r"([a-zA-Z])(\d+)(?=$|[/*.)])" => s"\1^\2")

    try
        u = uparse(s; unit_context=_CDS_UNIT_CONTEXT)
        return (; unit = u, valuefn)
    catch
        @warn "Could not parse CDS unit string" unit_string=s
        return (; unit = NoUnits, valuefn)
    end
end
