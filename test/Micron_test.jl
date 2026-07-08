# Tests that μm function as a unit

using PSSFSS
using Test


@testset "micron" begin
    @test 1u"mm" == 1000u"μm"
end
