# Tests that translating geometry does not break the simulation

using PSSFSS
using Test
using Unitful

# Test that the simulation runs normally for a variety of translations and mesh densities

function run_sim(ntri, shift)
    period = 20.0 # um
    l1 = 0.8*period
    l2 = 0.1*period

    flist = [5u"THz"]
    steering = (phi=0, theta=0)
    sheet = loadedcross(s1=[period,0.0], s2=[0.0,period], L1=l1, L2=l2, w=l2, units=μm, ntri=ntri, class='M', dx=shift*period, dy=shift*period)
    strata = [Layer(), sheet, Layer()]
    try
        results = analyze(strata, flist, steering, showprogress=false,resultfile=devnull, logfile=devnull)
        return true
    catch e 
        return false
    end
end

@testset "translate" for ntri in [100, 200, 400], shift in [0.0, 0.1, 0.25, 0.45, 0.50]
    @test run_sim(ntri, shift)
end