# MagnetoTransport.jl

[![CI](https://github.com/ChristopheBerthod/MagnetoTransport.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ChristopheBerthod/MagnetoTransport.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/ChristopheBerthod/MagnetoTransport.jl/graph/badge.svg?token=cXaZZi9hdM)](https://codecov.io/gh/ChristopheBerthod/MagnetoTransport.jl)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/ChristopheBerthod/MagnetoTransport.jl/blob/main/LICENSE)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://ChristopheBerthod.github.io/MagnetoTransport.jl/dev)

Solution of the electronic linear magneto-transport problem in two dimensions with a local energy-dependent self-energy, written in [Julia](https://julialang.org/). See the [Documentation](https://ChristopheBerthod.github.io/MagnetoTransport.jl/dev) for a definition of the problem.

## Dependencies

This package requires the  [Piecewise](https://github.com/ChristopheBerthod/Piecewise.jl) suite, which in turn makes use of special functions provided ny [HypergeometricFunctions.jl](https://github.com/JuliaMath/HypergeometricFunctions.jl), and [PolyLog](https://github.com/Expander/PolyLog.jl), as well as [LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl).

## Reference

G. Morpurgo, L. Rademaker, C. Berthod, and T. Giamarchi, [Physical Review Research **6**, 013112 (2024)](https://doi.org/10.1103/PhysRevResearch.6.013112).