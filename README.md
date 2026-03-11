# VOUnits.jl

Parse astronomical unit strings — as defined by the [IVOA VOUnit standard](https://www.ivoa.net/documents/VOUnits/) and related formats — into [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) units. Supported formats:

- **VOUnit** — IVOA standard (e.g. `m**-2`, `erg.s**-1.cm**-2`)
- **CDS** — used by VizieR and CDS catalogs (e.g. `km.s-1`, `[K]`, `10+3J/m/s/kpc2`)
- **FITS** — FITS file headers (similar to VOUnit)

This package is currently used in [VOTables.jl](https://github.com/JuliaAPlavin/VOTables.jl), [AstroTables.jl](https://github.com/JuliaAPlavin/AstroTables.jl), and [EnhancedCSV.jl](https://github.com/JuliaAPlavin/EnhancedCSV.jl) for units handling.

## Usage

By default, `parse_unit` accepts accepts all format variants and is liberal about which syntax you throw at it:

```julia
julia> using VOUnits

julia> parse_unit("km/s")
(unit = km s⁻¹, valuefn = identity)

julia> parse_unit("km.s-1")       # CDS-style dot multiplication and bare-sign powers
(unit = km s⁻¹, valuefn = identity)

julia> parse_unit("erg.s**-1.cm**-2")  # VOUnit/FITS ** power syntax
(unit = erg cm⁻² s⁻¹, valuefn = identity)

julia> parse_unit("solMass/yr")
(unit = M⊙ yr⁻¹, valuefn = identity)
```

Log/bracket notation (common in CDS catalogs for quantities in dex) returns an `exp10` transform:

```julia
julia> parse_unit("[K]")
(unit = K, valuefn = exp10)

julia> parse_unit("log(Hz)")
(unit = Hz, valuefn = exp10)
```

Convert back from Unitful units to a VOUnit-format string:

```julia
julia> unit_string(u"km/s")
"km.s**-1"

julia> unit_string(u"Msun")
"solMass"
```

If you need format-specific parsing, pass a format specifier: `parse_unit(s, CDS())`, `parse_unit(s, VOUnit())`, or `parse_unit(s, FITS())`. There are barely any differences now, but more specific checks can be added in the future.
