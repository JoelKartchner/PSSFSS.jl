using PSSFSS
using Test

sh1 = rectstrip(Lx=1, Ly=1.0, Nx=1, Ny=1, Px=1, Py=1, units=inch)

@testset "recstrip" begin
    @test length(sh1.ρ) == 4
    @test sh1.ρ[2] - sh1.ρ[1] == [1, 0]
    @test sh1.ρ[3] - sh1.ρ[1] == [0, 1]
    @test sh1.ρ[4] - sh1.ρ[1] == [1, 1]
    @test sh1.e1 == [1, 3, 1, 2, 1]
    @test sh1.e2 == [2, 4, 3, 4, 4]
    @test sh1.fe == reshape([2, 3, 5, 4, 5, 1], 3, 2)
    @test sh1.fv == reshape([1, 4, 3, 1, 2, 4], 3, 2)
    @test sh1.fufp
    @test sh1.dx == sh1.dy == 0
    @test sh1.rot == 0
    @test sh1.units == inch
    @test sh1.class == 'J'
    @test sh1.s₁ == [1, 0]
    @test sh1.s₂ == [0, 1]
    @test sh1.β₁ ≈ 2π .* [1, 0]
    @test sh1.β₂ ≈ 2π .* [0, 1]
    @test sh1.Zs == 0
end

Z = 0.2 + 0.3im
sh2 = rectstrip(Lx=1, Ly=1.0, Nx=1, Ny=1, Px=1, Py=1, units=inch, Zsheet=Z)
@testset "Zsheet" begin
    @test Z == sh2.Zs
end

R = 0.25
sh3 = rectstrip(Lx=1, Ly=1.0, Nx=1, Ny=1, Px=1, Py=1, units=inch, Rsheet=R)
@testset "Rsheet" begin
    @test R == sh3.Zs
end

@testset "badZsheet" begin
    @test_throws ErrorException rectstrip(Lx=1, Ly=1.0, Nx=1, Ny=1, Px=1, Py=1, units=inch, Zsheet=-0.2+0.2im)
end
