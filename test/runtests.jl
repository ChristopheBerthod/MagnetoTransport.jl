using MagnetoTransport
using Piecewise, PiecewiseHilbert, PiecewiseLorentz
using Test

@testset "number and chemical_potential" begin

    HoN0 = HilbertTransform(PiecewiseFunction(:even, Piece((0, 1), POLY, [1])))

    @test number(0, 0, HoN0) == 0.9999999982719904
    @test number(0, 0, HoN0, tol=1e-16) == 1.0000000000000029
    @test number(0, 1, HoN0, tol=1e-16) == 1.0000000000000029
    @test number(-10, 0, HoN0) == 0.0
    @test number(10, 0, HoN0) == 2.0
    @test number(0, 0, HoN0, ε -> -im) == 0.9999999977430993
    @test number(-10, 0, HoN0, ε -> -im) == 0.06366027892951089
    @test number(10, 0, HoN0, ε -> -im) == 1.9363397172066177
    @test number(0, 1, HoN0, ε -> -im) == 0.9999999977430993
    
    @test chemical_potential(1, 0, HoN0, tol=1e-4) == -1.698677792014653e-5
    @test chemical_potential(1, 0, HoN0, μ0=1, tol=1e-4) == -1.6986783310020074e-5
    @test chemical_potential(1, 0, HoN0, μ0=-1, tol=1e-4) == -1.6986783310020074e-5

end

@testset "σ₀, σ₁, and RH" begin

    HoN0 = HilbertTransform(PiecewiseFunction(:even, Piece((0, 1), POLY, [1])))
    Φ₀ = PiecewiseFunction(:even, Piece((0, 1), POLY, [1, -1]))
    Φ₁ = 2 * π^2 / 3 * PiecewiseFunction(:odd,
        [Piece((0, 0.5), POLY, [0, 1]), Piece((0.5, 1), POLY, [1, -1])])
    L2oPhi0 = LorentzTransform(Φ₀, 2)
    L3oPhi1 = LorentzTransform(Φ₁, 3)

    @test σ₀(1, 0, L2oPhi0, ε -> -im, HoN0) == 1.5707963267948963
    @test σ₀(1, 0, L2oPhi0, ε -> -im, HoN0, εmin=-1) == 1.5707963267948963
    @test σ₀(1, 0, L2oPhi0, ε -> -im, HoN0, εmax=1) == 1.5707963267948963
    @test σ₀(1, 0, L2oPhi0, ε -> -im, HoN0, μ0=-1) == 1.5707963267948963
    @test σ₀(1, 0, L2oPhi0, ε -> -im, HoN0, tol=1e-4) == 1.5707963267948963
    @test σ₀(1, 1, L2oPhi0, ε -> -im, HoN0) == 0.6692950293036449
    @test σ₀(1, 1, L2oPhi0, ε -> -im, HoN0, εmin=-1) == 0.616025642087738
    @test σ₀(1, 1, L2oPhi0, ε -> -im, HoN0, εmax=1) == 0.6160256416170858
    @test σ₀(1, 1, L2oPhi0, ε -> -im, HoN0, μ0=-1) == 0.669295029303645
    @test σ₀(1, 1, L2oPhi0, ε -> -im, HoN0, tol=1e-4) == 0.6692950293068127
    @test σ₀(1, 0, L2oPhi0, ε -> -im) == 0.6435011087932846
    @test σ₀(1, 0, L2oPhi0, ε -> -im, εmin=-1) == 0.6435011087932846
    @test σ₀(1, 0, L2oPhi0, ε -> -im, εmax=1) == 0.6435011087932846
    @test σ₀(1, 0, L2oPhi0, ε -> -im, tol=1e-4) == 0.6435011087932846
    @test σ₀(1, 1, L2oPhi0, ε -> -im) == 0.5657735906762745
    @test σ₀(1, 1, L2oPhi0, ε -> -im, εmin=-1) == 0.5398517772930664
    @test σ₀(1, 1, L2oPhi0, ε -> -im, εmax=1) == 0.4850732554158676
    @test σ₀(1, 1, L2oPhi0, ε -> -im, tol=1e-4) == 0.5657735917212093

    @test σ₁(1, 0, L3oPhi1, ε -> -im, HoN0) == 1.990652773064181e-8
    @test σ₁(1, 0, L3oPhi1, ε -> -im, HoN0, εmin=-1) == 1.990652773064181e-8
    @test σ₁(1, 0, L3oPhi1, ε -> -im, HoN0, εmax=1) == 1.990652773064181e-8
    @test σ₁(1, 0, L3oPhi1, ε -> -im, HoN0, μ0=-1) == 1.9906522936752167e-8
    @test σ₁(1, 0, L3oPhi1, ε -> -im, HoN0, tol=1e-4) == 8.012683160043856e-9
    @test σ₁(1, 1, L3oPhi1, ε -> -im, HoN0) == 1.6324236969273102e-9
    @test σ₁(1, 1, L3oPhi1, ε -> -im, HoN0, εmin=-1) == 0.06682971964324791
    @test σ₁(1, 1, L3oPhi1, ε -> -im, HoN0, εmax=1) == -0.06682971693862014
    @test σ₁(1, 1, L3oPhi1, ε -> -im, HoN0, μ0=-1) == 1.632381551962753e-9
    @test σ₁(1, 1, L3oPhi1, ε -> -im, HoN0, tol=1e-4) == 6.570754696208626e-10
    @test σ₁(1, 0, L3oPhi1, ε -> -im) == 1.0429289227284044
    @test σ₁(1, 0, L3oPhi1, ε -> -im, εmin=-1) == 1.0429289227284044
    @test σ₁(1, 0, L3oPhi1, ε -> -im, εmax=1) == 1.0429289227284044
    @test σ₁(1, 0, L3oPhi1, ε -> -im, tol=1e-4) == 1.0429289227284044
    @test σ₁(1, 1, L3oPhi1, ε -> -im) == 0.17322935899094005
    @test σ₁(1, 1, L3oPhi1, ε -> -im, εmin=-1) == 0.20642633236588082
    @test σ₁(1, 1, L3oPhi1, ε -> -im, εmax=1) == 0.07640628017917533
    @test σ₁(1, 1, L3oPhi1, ε -> -im, tol=1e-4) == 0.17322935898931469

    RH(1, 0, L2oPhi0, L3oPhi1, ε -> -im, HoN0) == 8.067811807510598e-9
    RH(1, 0, L2oPhi0, LorentzTransform(PiecewiseFunction(:even,
        Piece((0, 1), POLY, [2, -2])), 2), L3oPhi1, ε -> -im, HoN0) == 4.033905903755299e-9
    RH(1, 0, L2oPhi0, L3oPhi1, ε -> -im) == 2.518582100162326
    RH(1, 0, L2oPhi0, LorentzTransform(PiecewiseFunction(:even,
        Piece((0, 1), POLY, [2, -2])), 2), L3oPhi1, ε -> -im) == 1.259291050081163

end
