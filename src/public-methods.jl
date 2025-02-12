
"""
    number(μ, T, HoN0, Σ [, tol])

Particle number for chemical potential `μ` and temperature `T` (energy units).
`HoN0` is a [`HilbertTransform`](@ref) object holding the Hilbert transform of
the density of states. `Σ(::Real)` is a function giving the self-energy. `tol`
sets the relative tolerance of the quadratures (default `1e-8`).
"""
function number(μ::Real, T::Real, HoN0::HilbertTransform, Σ::Function; tol=1e-8)
    # Interacting DOS
    N(E) = -imag(HoN0(μ + E - Σ(μ + E))) / π
    # Zero-temperature contribution
    lb, ub = support(HoN0.f)
    if N(0.0) < 100 * eps(Float64)
        # Case where the chemical potential is outside the band
        n = μ < lb ? 0.0 : HoN0.moments[1]
    else
        # Find lower bound of integral
        if -imag(Σ(-Inf)) > eps(Float64)
            Eₘᵢₙ = -Inf
        else
            Eₘᵢₙ = lb - μ
            while N(Eₘᵢₙ) > 100 * eps(Float64)
                Eₘᵢₙ -= 0.1 * (ub - lb)
            end
        end
        # We keep the error tested by rtol above 100 * eps(Float64).
        # To this end, we estimate the integral with rtol=1e-2...
        int = quadgk(E -> N(E), Eₘᵢₙ, 0, maxevals=1e3,
            rtol=1e-2, atol=eps(Float64))
        # ... and we then enforce int[1] * rtol > 100 * eps(Float64)
        int = quadgk(E -> N(E), Eₘᵢₙ, 0, maxevals=1e4,
            rtol=max(tol, 100 * eps(Float64) / int[1]), atol=eps(Float64))
        # Warnings are discarded for values of order 10^-8 or less
        if abs(int[1]) > eps(Float32) && int[2] / int[1] > tol
            e = @sprintf("%.1e", int[2] / int[1])
            @warn "Integral 1 not converged in number($(μ), $(T), HoN0, Σ; "*
                "tol=$(tol)) [error = $(e)]"
        end
        n = int[1]
    end
    # Finite-temperature correction
    if T > 0
        # like before...
        int = quadgk(u -> T * (N(T * u) - N(-T * u)) * _Fermi(u), 0, 10, 30, 100,
            rtol=1e-2, atol=eps(Float64))
        # ...
        int = quadgk(u -> T * (N(T * u) - N(-T * u)) * _Fermi(u), 0, 10, 30, 100,
            rtol=max(tol, 100 * eps(Float64) / int[1]), atol=eps(Float64))
        # Warnings are discarded for values of order 10^-8 or less
        if abs(int[1]) > eps(Float32) && int[2] / int[1] > tol
            e = @sprintf("%.1e", int[2] / int[1])
            @warn "Integral 2 not converged in number($(μ), $(T), HoN0, Σ; "*
                "tol=$(tol)) [error = $(e)]"
        end
        n += int[1]
    end
    return n
end


"""
    number(μ, T, HoN0 [, tol])

Particle number for chemical potential `μ` and temperature `T` (energy units) in
the noninteracting case (``\\Sigma(\\varepsilon)=-i0^+``). `HoN0` is a
[`HilbertTransform`](@ref) object holding the Hilbert transform of the density
of states. `tol` sets the relative tolerance of the quadratures (default
`1e-8`).
"""
number(μ::Real, T::Real, HoN0::HilbertTransform; tol=1e-8) =
    number(μ, T, HoN0, E -> -im * eps(Float64), tol=tol)


"""
    chemical_potential(n, T, HoN0, Σ [, μ0] [, tol])

Chemical potential for density `n` and temperature `T` (energy units). `HoN0` is
a [`HilbertTransform`](@ref) object holding the Hilbert transform of the density
of states. `Σ(::Real)` is a function giving the self-energy. `μ0` is an
optional initial estimate. `tol` sets the relative tolerance of the quadratures
and root finding (default `1e-8`).
"""
function chemical_potential(n::Real, T::Real, HoN0::HilbertTransform,
    Σ::Function; μ0=missing, tol=1e-8)
    # Initial estimate
    lb, ub = support(HoN0.f)
    μ₀ = ismissing(μ0) ? lb + n * (ub - lb) / HoN0.moments[1] : μ0
    # DOS at μ₀
    N0 = -imag(HoN0(μ₀ - Σ(μ₀))) / π
    # Initial width of the root-searching interval: 10% change in n
    w = 0.1 * n / (N0 > 100 * eps(Float64) ? N0 : HoN0.moments[1] / (ub - lb))
    # Root finding
    zero(μ) = number(μ, T, HoN0, Σ, tol=tol) - n
    z0 = zero(μ₀) 
    if z0 > 0
        z1 = zero(μ₀ - w)
        while z1 > 0
            w *= 10
            z1 = zero(μ₀ - w)
        end
        return _find_zero_Ridders(zero, (μ₀ - w, μ₀); values=(z1, z0), tol=tol)
    else
        z1 = zero(μ₀ + w)
        while z1 < 0
            w *= 10
            z1 = zero(μ₀ + w)
        end
        return _find_zero_Ridders(zero, (μ₀, μ₀ + w); values=(z0, z1), tol=tol)
    end
end


"""
    chemical_potential(n, T, HoN0 [, μ0] [, tol])

Chemical potential for density `n` and temperature `T` (energy units) in the
noninteracting case (``\\Sigma(\\varepsilon)=-i0^+``). `HoN0` is a
[`HilbertTransform`](@ref) object holding the Hilbert transform of the density
of states. `μ0` is an optional initial estimate. `tol` sets the relative
tolerance of the quadratures and root finding (default `1e-8`).
"""
chemical_potential(n::Real, T::Real, HoN0::HilbertTransform; μ0=missing, tol=1e-8) =
    chemical_potential(n, T, HoN0, E -> -im * eps(Float64); μ0=μ0, tol=tol)


"""
    σ₀(n, T, L2oPhi0, Σ, HoN0 [, εmin] [, εmax] [, μ0] [, tol])

Normal conductivity in units of ``e^2/h`` for density `n` and temperature `T`
(energy units). `L2oPhi0` is a [`LorentzTransform`](@ref) object holding the
``L^2`` transform of the transport function divided by ``(e/\\hbar)^2``.
`Σ(::Real)` is a function giving the self-energy. `HoN0` is a
[`HilbertTransform`](@ref) object holding the Hilbert transform of the density
of states. `εmin` and `εmax` are optional limits of integration at finite-``T``
(default `-Inf` and `+Inf`). `μ0` is an optional initial estimate for the
chemical potential. `tol` sets the relative tolerance of the quadratures and
root finding (default `1e-8`).
"""
function σ₀(n::Real, T::Real, L2oPhi0::LorentzTransform, Σ::Function,
    HoN0::HilbertTransform; εmin=missing, εmax=missing, μ0=missing, tol=1e-8)
    L2oPhi0.m == 2 || throw(ArgumentError(
        "The argument passed as L2oPhi0 is not an order-2 Lorentz transform."))
    μ = chemical_potential(n, T, HoN0, Σ; μ0=μ0, tol=tol)
    if T ≈ 0
        return L2oPhi0(μ, Σ(μ)) * 2 * π^2
    else
        bounds = μ .+ [-100, -30, -10, 0, 10, 30, 100] .* T
        !ismissing(εmin) && (bounds = vcat(εmin, filter(x -> x > εmin, bounds)))
        !ismissing(εmax) && (bounds = vcat(filter(x -> x < εmax, bounds), εmax))
        bounds .= (bounds .- μ) ./ T
        int = quadgk(u -> L2oPhi0(T * u + μ, Σ(T * u + μ)) * _dFermi(u),
            bounds..., rtol=tol)
        if int[2] / int[1] > tol
            e = @sprintf("%.1e", int[2] / int[1])
            @warn "Integral not converged in σ₀($(n), $(T), L2oPhi0, Σ, HoN0; " *
                "εmin=$(εmin), μ0=$(μ0), tol=$(tol)) [error = $(e)]"
        end
        return int[1]  * 2 * π^2
    end
end


"""
    σ₀(μ, T, L2oPhi0, Σ [, εmin] [, εmax] [, tol])

Normal conductivity in units of ``e^2/h`` for chemical potential `μ` and
temperature `T` (energy units). `L2oPhi0` is a [`LorentzTransform`](@ref) object
holding the ``L^2`` transform of the transport function divided by
``(e/\\hbar)^2``. `Σ(::Real)` is a function giving the self-energy. `εmin` and
`εmax` are optional limits of integration at finite-``T`` (default `-Inf` and
`+Inf`). `tol` sets the relative tolerance of the quadrature (default `1e-8`). 
"""
function σ₀(μ::Real, T::Real, L2oPhi0::LorentzTransform, Σ::Function;
    εmin=missing, εmax=missing, tol=1e-8)
    L2oPhi0.m == 2 || throw(ArgumentError(
        "The argument passed as L2oPhi0 is not an order-2 Lorentz transform."))
    if T ≈ 0
        return L2oPhi0(μ, Σ(μ)) * 2 * π^2
    else
        bounds = μ .+ [-100, -30, -10, 0, 10, 30, 100] .* T
        !ismissing(εmin) && (bounds = vcat(εmin, filter(x -> x > εmin, bounds)))
        !ismissing(εmax) && (bounds = vcat(filter(x -> x < εmax, bounds), εmax))
        bounds .= (bounds .- μ) ./ T
        int = quadgk(u -> L2oPhi0(T * u + μ, Σ(T * u + μ)) * _dFermi(u),
            bounds..., rtol=tol)
        if int[2] / int[1] > tol
            e = @sprintf("%.1e", int[2] / int[1])
            @warn "Integral not converged in σ₀($(μ), $(T), L2oPhi0, Σ; " *
                "εmin=$(εmin), tol=$(tol)) [error = $(e)]"
        end
        return int[1]  * 2 * π^2
    end
end


"""
    σ₁(n, T, L3oPhi1, Σ, HoN0 [, εmin] [, εmax] [, μ0] [, tol])

Hall conductivity divided by ``B|e|^3/h^2`` for density `n` and temperature `T`
(energy units). `L3oPhi1` is a [`LorentzTransform`](@ref) object holding the
``L^3`` transform of the transport function divided by ``(|e|/\\hbar)^3``.
`Σ(::Real)` is a function giving the self-energy. `HoN0` is a
[`HilbertTransform`](@ref) object holding the Hilbert transform of the density
of states. `εmin` and `εmax` are optional limits of integration at finite-``T``
(default `-Inf` and `+Inf`). `μ0` is an optional initial estimate for the
chemical potential. `tol` sets the relative tolerance of the quadratures and
root finding (default `1e-8`).
"""
function σ₁(n::Real, T::Real, L3oPhi1::LorentzTransform, Σ::Function,
    HoN0::HilbertTransform; εmin=missing, εmax=missing, μ0=missing, tol=1e-8)
    L3oPhi1.m == 3 || throw(ArgumentError(
        "The argument passed as L3oPhi1 is not an order-3 Lorentz transform."))
    μ = chemical_potential(n, T, HoN0, Σ; μ0=μ0, tol=tol)
    if T ≈ 0
        return L3oPhi1(μ, Σ(μ)) * 4 * π^2
    else
        bounds = μ .+ [-100, -30, -10, 0, 10, 30, 100] .* T
        !ismissing(εmin) && (bounds = vcat(εmin, filter(x -> x > εmin, bounds)))
        !ismissing(εmax) && (bounds = vcat(filter(x -> x < εmax, bounds), εmax))
        bounds .= (bounds .- μ) ./ T
        int = quadgk(u -> L3oPhi1(T * u + μ, Σ(T * u + μ)) * _dFermi(u),
            bounds..., rtol=tol)
        if int[2] / int[1] > tol
            e = @sprintf("%.1e", int[2] / int[1])
            @warn "Integral not converged in σ₁($(n), $(T), L3oPhi1, Σ, HoN0; " *
                "εmin=$(εmin), μ0=$(μ0), tol=$(tol)) [error = $(e)]"
        end
        return int[1]  * 4 * π^2
    end
end


"""
    σ₁(μ, T, L3oPhi1, Σ [, εmin] [, εmax] [, tol])

Hall conductivity divided by ``B|e|^3/h^2`` for chemical potential `μ` and
temperature `T` (energy units). `L3oPhi1` is a [`LorentzTransform`](@ref) object
holding the ``L^3`` transform of the transport function divided by
``(|e|/\\hbar)^3``. `Σ(::Real)` is a function giving the self-energy. `εmin` and
`εmax` are optional limits of integration at finite-``T`` (default `-Inf` and
`+Inf`). `tol` sets the relative tolerance of the quadrature (default `1e-8`).
"""
function σ₁(μ::Real, T::Real, L3oPhi1::LorentzTransform, Σ::Function;
    εmin=missing, εmax=missing, tol=1e-8)
    L3oPhi1.m == 3 || throw(ArgumentError(
        "The argument passed as L3oPhi1 is not an order-3 Lorentz transform."))
    if T ≈ 0
        return L3oPhi1(μ, Σ(μ)) * 4 * π^2
    else
        bounds = μ .+ [-100, -30, -10, 0, 10, 30, 100] .* T
        !ismissing(εmin) && (bounds = vcat(εmin, filter(x -> x > εmin, bounds)))
        !ismissing(εmax) && (bounds = vcat(filter(x -> x < εmax, bounds), εmax))
        bounds .= (bounds .- μ) ./ T
        int = quadgk(u -> L3oPhi1(T * u + μ, Σ(T * u + μ)) * _dFermi(u),
            bounds..., rtol=tol)
        if int[2] / int[1] > tol
            e = @sprintf("%.1e", int[2] / int[1])
            @warn "Integral not converged in σ₁($(μ), $(T), L3oPhi1, Σ; " *
                "εmin=$(εmin), tol=$(tol)) [error = $(e)]"
        end
        return int[1]  * 4 * π^2
    end
end


"""
    RH(n, T, L2oPhi01 [, L2oPhi02], L3oPhi1, Σ, HoN0 [, εmin] [, εmax] [, μ0] [, tol])

Hall constant multiplied by ``|e|`` for density `n` and temperature `T` (energy
units). `L2oPhi01`, `L2oPhi02`, and `L3oPhi1` are [`LorentzTransform`](@ref)
objects holding the ``L^2`` and ``L^3`` transforms of the transport functions.
`Σ(::Real)` is a function giving the self-energy. `HoN0` is a
[`HilbertTransform`](@ref) object holding the Hilbert transform of the density
of states. `εmin` and `εmax` are optional limits of integration at finite-``T``
(default `-Inf` and `+Inf`). `μ0` is an optional initial estimate for the
chemical potential. `tol` sets the relative tolerance of the quadratures and
root finding (default `1e-8`).
"""
function RH(n::Real, T::Real, L2oPhi01::LorentzTransform,
    L2oPhi02::LorentzTransform, L3oPhi1::LorentzTransform, Σ::Function,
    HoN0::HilbertTransform; εmin=missing, εmax=missing, μ0=missing, tol=1e-8)
    μ = chemical_potential(n, T, HoN0, Σ; μ0=μ0, tol=tol)
    if L2oPhi01 === L2oPhi02
        return σ₁(μ, T, L3oPhi1, Σ; εmin=εmin, εmax=εmax, tol=tol) /
            σ₀(μ, T, L2oPhi01, Σ; εmin=εmin, εmax=εmax, tol=tol)^2
    else
        return σ₁(μ, T, L3oPhi1, Σ; εmin=εmin, εmax=εmax, tol=tol) / (
            σ₀(μ, T, L2oPhi01, Σ; εmin=εmin, εmax=εmax, tol=tol) *
            σ₀(μ, T, L2oPhi02, Σ; εmin=εmin, εmax=εmax, tol=tol))
    end
end

RH(n::Real, T::Real, L2oPhi0::LorentzTransform, L3oPhi1::LorentzTransform,
    Σ::Function, HoN0::HilbertTransform;
    εmin=missing, εmax=missing, μ0=missing, tol=1e-8) =
    RH(n, T, L2oPhi0, L2oPhi0, L3oPhi1, Σ, HoN0; εmin=εmin, εmax=εmax, μ0=μ0, tol=tol)


"""
    RH(μ, T, L2oPhi01 [, L2oPhi02], L3oPhi1, Σ [, εmin] [, εmax] [, tol])

Hall constant multiplied by ``|e|`` for chemical potential `μ` and temperature
`T` (energy units). `L2oPhi01`, `L2oPhi02`, and `L3oPhi1` are
[`LorentzTransform`](@ref) objects holding the ``L^2`` and ``L^3`` transforms of
the transport functions. `Σ(::Real)` is a function giving the self-energy.
`εmin` and `εmax` are optional limits of integration at finite-``T`` (default
`-Inf` and `+Inf`). `tol` sets the relative tolerance of the quadratures
(default `1e-8`).
"""
function RH(μ::Real, T::Real, L2oPhi01::LorentzTransform,
    L2oPhi02::LorentzTransform, L3oPhi1::LorentzTransform, Σ::Function;
    εmin=missing, εmax=missing, tol=1e-8)
    if L2oPhi01 === L2oPhi02
        return σ₁(μ, T, L3oPhi1, Σ; εmin=εmin, εmax=εmax, tol=tol) /
        σ₀(μ, T, L2oPhi01, Σ; εmin=εmin, εmax=εmax, tol=tol)^2
    else
        return σ₁(μ, T, L3oPhi1, Σ; εmin=εmin, εmax=εmax, tol=tol) / (
            σ₀(μ, T, L2oPhi01, Σ; εmin=εmin, εmax=εmax, tol=tol) *
            σ₀(μ, T, L2oPhi02, Σ; εmin=εmin, εmax=εmax, tol=tol))
    end
end

RH(μ::Real, T::Real, L2oPhi0::LorentzTransform, L3oPhi1::LorentzTransform,
    Σ::Function; εmin=missing, εmax=missing, tol=1e-8) =
    RH(μ, T, L2oPhi0, L2oPhi0, L3oPhi1, Σ, εmin=εmin, εmax=εmax, tol=tol)
