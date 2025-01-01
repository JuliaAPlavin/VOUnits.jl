using TestItems
using TestItemRunner
@run_package_tests


@testitem "_" begin
    import Aqua
    Aqua.test_all(AstroASCIITables; ambiguities=false)
    Aqua.test_ambiguities(AstroASCIITables)

    import CompatHelperLocal as CHL
    CHL.@check()
end
