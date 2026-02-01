# Tests for CDS unit string → Unitful conversion.
# Ported from astropy: astropy/units/tests/test_format.py
# test_cds_grammar (lines 91-138), test_cds_grammar_fail (lines 141-170),
# test_cds_dimensionless, test_cds_log10_dimensionless, test_percent,
# test_scaled_dimensionless.

@testitem "cds_to_unitful grammar" begin
    using Unitful
    using UnitfulAstro
    using UnitfulAngles

    _unit(s) = cds_to_unitful(s).unit
    _valuefn(s) = cds_to_unitful(s).valuefn

    # ── test_cds_grammar: successful parses ──────────────────────────────
    # Each line corresponds to an entry in astropy's test_cds_grammar parametrize.

    @test _unit("0.1nm") == 0.1u"nm"                          # (["0.1nm"], u.AA)
    @test _unit("mW/m2") == u"mW/m^2"                         # (["mW/m2"], u.Unit(u.erg/u.cm**2/u.s))
    @test _unit("mW/(m2)") == u"mW/m^2"                       # (["mW/(m2)"], ...)
    @test _unit("km/s") == u"km/s"                             # (["km/s","km.s-1"], u.km/u.s)
    @test _unit("km.s-1") == u"km/s"                           # same pair
    @test _unit("km/s/Mpc") == u"km/s/Mpc"                    # (["km/s/Mpc"], ...)
    @test _unit("km/(s.Mpc)") == u"km/s/Mpc"                  # (["km/(s.Mpc)"], ...)
    @test _unit("10+3J/m/s/kpc2") == 1e3u"J/m/s/kpc^2"       # (["10+3J/m/s/kpc2"], u.Unit(1e3*u.W/(u.m*u.kpc**2)))
    @test _unit("10pix/nm") == 10u"nm^-1"                       # pix removed as dimensionless
    @test _unit("1.5x10+11m") == 1.5e11u"m"                   # (["1.5x10+11m"], ...)
    @test _unit("1.5×10+11/m") == 1.5e11u"m^-1"               # (["1.5×10+11/m"], ...)
    @test _unit("/s") == u"s^-1"                               # (["/s"], u.s**-1)
    @test _unit("m2") == u"m^2"                                # (["m2"], u.m**2)
    @test _unit("10+21m") == 1e21u"m"                          # (["10+21m"], ...)
    @test _unit("2.54cm") == 2.54u"cm"                         # (["2.54cm"], ...)
    @test _unit("20%") == 20u"percent"                         # (["20%"], 0.20*u.dimensionless_unscaled)
    @test _unit("10+9") == 1e9 * NoUnits                       # (["10+9"], 1.0e9*u.dimensionless_unscaled)
    @test _unit("2x10-9") == 2e-9 * NoUnits                    # (["2x10-9"], ...)
    @test _unit("---") == NoUnits                              # (["---"], u.dimensionless_unscaled)
    @test_broken _unit("ma") != NoUnits                        # "a" (year) has no SI-prefixed forms in Unitful
    @test_broken _unit("mAU") != NoUnits                       # AU has no SI-prefixed forms in Unitful
    @test_broken _unit("uarcmin") != NoUnits                   # micro-arcminute not in UnitfulAngles
    @test _unit("uarcsec") == u"μas"                           # micro-arcsecond = μas in UnitfulAngles
    @test_broken _unit("kbarn") != NoUnits                     # barn not in Unitful
    @test_broken _unit("Gbit") != NoUnits                      # bit not in Unitful
    @test_broken _unit("Gibit") != NoUnits                     # binary prefix not in Unitful
    @test_broken _unit("kbyte") != NoUnits                     # byte not in Unitful
    @test_broken _unit("mRy") != NoUnits                       # Rydberg not in Unitful
    @test_broken _unit("mmag") != NoUnits                      # mag has no SI-prefixed forms in Unitful
    @test _unit("Mpc") == u"Mpc"                               # (["Mpc"], u.Mpc)
    @test _unit("Gyr") == u"Gyr"                               # (["Gyr"], u.Gyr)
    @test _unit("°") == u"°"                                   # (["°"], u.degree)
    @test _unit("°/s") == u"°/s"                               # (["°/s"], u.degree/u.s)
    @test _unit("Å") == u"Å"                                   # (["Å"], u.AA)
    @test _unit("Å/s") == u"Å/s"                               # (["Å/s"], u.AA/u.s)
    @test _unit("\\h") == Unitful.h                                    # Planck constant via CDS \\h notation
    @test _unit("eps0/mu0") != NoUnits                          # physical constants via name mapping
    @test _unit("a0.s") == NoUnits                              # Bohr radius not in Unitful; mapped to bohrRadius → warns
end

@testitem "cds_to_unitful dex" begin
    using Unitful
    using UnitfulAstro
    using UnitfulAngles

    _unit(s) = cds_to_unitful(s).unit
    _valuefn(s) = cds_to_unitful(s).valuefn

    # ── test_cds_grammar: bracket/dex notation ────────────────────────────
    # CDS [unit] = dex (log10) of unit. We return the base unit + exp10 transform.
    @test _unit("[cm/s2]") == u"cm/s^2"                        # (["[cm/s2]"], dex(u.cm/u.s**2))
    @test _valuefn("[cm/s2]") === exp10
    @test _unit("[K]") == u"K"                                 # (["[K]"], dex(u.K))
    @test _valuefn("[K]") === exp10
    @test _unit("[-]") == NoUnits                              # (["[-]"], dex(u.dimensionless_unscaled))
    @test _valuefn("[-]") === exp10

    # ── test_cds_dimensionless ────────────────────────────────────────────
    @test _unit("---") == NoUnits
    @test _valuefn("---") === identity

    # ── test_cds_log10_dimensionless ──────────────────────────────────────
    @test _unit("[-]") == NoUnits
    @test _valuefn("[-]") === exp10
end

@testitem "cds_to_unitful fail" begin
    using Unitful

    _unit(s) = cds_to_unitful(s).unit

    # ── test_cds_grammar_fail: these are invalid CDS unit strings ────────
    # In astropy, all of these raise ValueError.
    # Our function returns NoUnits + @warn instead of throwing.
    @test _unit("0.1 nm") == NoUnits       # space in unit
    @test _unit("solMass(3/2)") == NoUnits  # exponent in parens
    @test_broken _unit("km / s") == NoUnits  # spaces: Unitful tolerates spaces, astropy rejects
    @test _unit("km s-1") == NoUnits        # space instead of dot
    @test_broken _unit("km/s.Mpc-1") == NoUnits  # mixed dot and slash: our transform makes it valid
    @test_broken _unit("/s.Mpc") == NoUnits       # leading slash + dot: our transform makes it valid
    @test_broken _unit("pix0.1nm") == NoUnits      # pix removed → "0.1nm" parses; invalid CDS but unlikely
    @test _unit("pix/(0.1nm)") == 10u"nm^-1"        # pix removed, parses as /(0.1nm)
    @test_broken _unit("km*s") == NoUnits   # * not valid in CDS, but Unitful accepts it
    @test _unit("km**2") == NoUnits         # ** not valid in CDS
    @test _unit("5x8+3m") == NoUnits        # non-base-10 factor
    @test _unit("0.1---") == NoUnits        # prefix before ---
    @test _unit("---m") == NoUnits          # --- before unit
    @test _unit("m---") == NoUnits          # unit before ---
    @test _unit("--") == NoUnits            # not a valid delimiter
    @test _unit("0.1-") == NoUnits          # prefix before -
    @test _unit("-m") == NoUnits            # - before unit
    @test _unit("m-") == NoUnits            # trailing minus
    @test _unit("mag(s-1)") == NoUnits      # function-call syntax
    @test _unit("dB(mW)") == NoUnits        # dB notation
    @test _unit("dex(cm s-2)") == NoUnits   # explicit dex() syntax
    @test _unit("[--]") == NoUnits          # double dash in brackets
end

@testitem "cds_to_unitful name mappings" begin
    using Unitful
    using UnitfulAstro
    using UnitfulAngles

    _unit(s) = cds_to_unitful(s).unit

    # ── CDS → Unitful name mappings ──────────────────────────────────────
    @test _unit("h") == u"hr"              # CDS h = hour, NOT Planck constant
    @test _unit("min") == u"minute"
    @test _unit("sec") == u"s"
    @test _unit("arcmin") == u"arcminute"
    @test _unit("arcsec") == u"arcsecond"
    @test _unit("arcm") == u"arcminute"
    @test _unit("arcs") == u"arcsecond"
    @test _unit("solMass") == u"Msun"
    @test _unit("solRad") == u"Rsun"
    @test _unit("solLum") == u"Lsun"
    @test _unit("Msun") == u"Msun"
    @test _unit("Rsun") == u"Rsun"
    @test _unit("Lsun") == u"Lsun"
    @test _unit("jovMass") == u"Mjup"
    @test _unit("Mjup") == u"Mjup"
    @test _unit("Rjup") == u"Rjup"
    @test _unit("Ohm") == u"Ω"
    @test _unit("AA") == u"Å"
    @test _unit("Angstrom") == u"Å"
    @test _unit("Angstroem") == u"Å"
    @test _unit("lyr") == u"ly"
    @test _unit("marcsec") == u"mas"
    @test _unit("uarcsec") == u"μas"
    @test _unit("gauss") == u"Gauss"
    @test _unit("Mgeo") == u"Mearth"
    @test _unit("Rgeo") == u"Rearth"
    @test _unit("au") == u"AU"              # CDS accepts both AU and au
end

@testitem "cds_to_unitful extra" begin
    using Unitful
    using UnitfulAstro
    using UnitfulAngles

    _unit(s) = cds_to_unitful(s).unit
    _valuefn(s) = cds_to_unitful(s).valuefn

    # ── Simple units (direct Unitful support) ─────────────────────────────
    @test _unit("K") == u"K"
    @test _unit("eV") == u"eV"
    @test _unit("deg") == u"°"
    @test _unit("mag") == u"mag"
    @test _unit("yr") == u"yr"
    @test _unit("d") == u"d"
    @test _unit("s") == u"s"
    @test _unit("rad") == u"rad"
    @test _unit("sr") == u"sr"
    @test _unit("erg") == u"erg"
    @test _unit("Jy") == u"Jy"
    @test _unit("pc") == u"pc"
    @test _unit("AU") == u"AU"
    @test _unit("bar") == u"bar"

    # ── SI prefixes ───────────────────────────────────────────────────────
    @test _unit("mJy") == u"mJy"
    @test _unit("kpc") == u"kpc"
    @test _unit("Mpc") == u"Mpc"
    @test _unit("Gyr") == u"Gyr"
    @test _unit("km") == u"km"

    # ── Decimal scale prefix ──────────────────────────────────────────────
    @test _unit("0.1pm") == 0.1u"pm"

    # ── Composite units from our test data ────────────────────────────────
    @test _unit("mW.m-2") == u"mW/m^2"
    @test _unit("m.s-2") == u"m/s^2"
    @test _unit("Msun/yr") == u"Msun/yr"
    @test _unit("Jy/sr") == u"Jy/sr"
    @test _unit("pc2") == u"pc^2"
    @test _unit("10+20cm-2") == 1e20u"cm^-2"

    # ── test_percent ──────────────────────────────────────────────────────
    @test _unit("%") == u"percent"

    # ── test_scaled_dimensionless ─────────────────────────────────────────
    @test _unit("10-4") == 1e-4 * NoUnits
    @test _unit("10+9") == 1e9 * NoUnits

    # ── Non-dex units have identity valuefn ───────────────────────────────
    @test _valuefn("K") === identity
    @test _valuefn("---") === identity
    @test _valuefn("km/s") === identity

    # ── Known dimensionless CDS units → replaced by 1 ──────────────────────
    @test _unit("ct") == 1 * NoUnits
    @test _unit("pix") == 1 * NoUnits
    @test _unit("ct/s") == 1u"s^-1"
    @test _unit("pix/nm") == 1u"nm^-1"
    @test _valuefn("ct") === identity

    # ── Dimensionless units with powers are replaced cleanly ─────────────
    @test _unit("Jy.beam-1") == 1u"Jy"       # beam^-1 is dimensionless
    @test _unit("Jy/beam") == 1.0u"Jy"       # Jy per beam → 1.0 Jy
    @test _unit("beam-1.Jy") == 1u"Jy"       # beam^-1 × Jy → 1 Jy
    @test _unit("Jy.beam-1.pix-1") == 1u"Jy" # multiple dimensionless with powers
    @test _unit("ct2/s") == 1u"s^-1"          # count^2 per second → 1 s^-1
    @test _unit("beam-1") == 1 * NoUnits      # bare dimensionless with power
end
