module MagnetoTransport

using QuadGK: quadgk
using Printf: @sprintf
using Piecewise, PiecewiseHilbert, PiecewiseLorentz

export

    # Methods
    number,
    chemical_potential,
    σ₀,
    σ₁,
    RH

include("private-methods.jl")
include("public-methods.jl")

end
