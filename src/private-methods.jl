
# Ridders method with bisection step
function _find_zero_Ridders(f::Function, bounds::Tuple{Real,Real}; 
    values::Tuple{Real,Real}=missing, tol::Real=1e-8)
    x1, x2 = bounds
    f1, f2 = ismissing(values) ? (f(x1), f(x2)) : values
    x3 = 0.0
    f1 * f2 < 0 || throw(ArgumentError("Not a bracketing interval."))
    while abs(x2 - x1) / (1 + abs(x3)) > tol
        x4 = (x1 + x2) / 2
        f4 = f(x4)
        x3 = x4 + (x4 - x1) * sign(f1) * f4 / sqrt(f4^2 - f1 * f2)
        f3 = f(x3); abs(f3) ≈ 0 && break
        if f1 * f3 < 0
            x2, f2 = x3, f3
        else
            x1, f1 = x3, f3
        end
        # Bisection step (to avoid possible cycles)
        x3 = (x1 + x2) / 2
        f3 = f(x3); abs(f3) ≈ 0 && break
        if f1 * f3 < 0
            x2, f2 = x3, f3
        else
            x1, f1 = x3, f3
        end
    end
    return x3
end


# Fermi distribution (without NaN risks)
_Fermi(u) = u < -100 ? 1.0 : u > 100 ? 0.0 : 1 / (exp(u) + 1)


# Derivative of the Fermi distribution (without NaN risks)
_dFermi(u) = abs(u) > 100 ? 0.0 : exp(u) / (exp(u) + 1)^2
