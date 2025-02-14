# MagnetoTransport.jl

[![CI](https://github.com/ChristopheBerthod/MagnetoTransport.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ChristopheBerthod/MagnetoTransport.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/ChristopheBerthod/MagnetoTransport.jl/graph/badge.svg?token=cXaZZi9hdM)](https://codecov.io/gh/ChristopheBerthod/MagnetoTransport.jl)
<!--[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/ChristopheBerthod/MagnetoTransport.jl/blob/main/LICENSE)-->
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://ChristopheBerthod.github.io/MagnetoTransport.jl/dev)

Solution of the electronic linear magneto-transport quantum problem in two dimensions with a local energy-dependent self-energy, written in [Julia](https://julialang.org/). See the [documentation](https://ChristopheBerthod.github.io/MagnetoTransport.jl/dev) for a precise definition of the problem.

### Dependencies

- [`Piecewise`, `PiecewiseHilbert`, `PiecewiseLorentz`](https://github.com/ChristopheBerthod/Piecewise.jl)

- [`QuadGK`](https://github.com/JuliaMath/QuadGK.jl)

### Installation

```julia
using Pkg
Pkg.add(url="https://github.com/ChristopheBerthod/MagnetoTransport.jl")
```

### Example

See the [example](https://ChristopheBerthod.github.io/MagnetoTransport.jl/dev/index.htm#Example) in the documentation.

### Reference

G. Morpurgo, L. Rademaker, C. Berthod, and T. Giamarchi, [Physical Review Research **6**, 013112 (2024)](https://doi.org/10.1103/PhysRevResearch.6.013112).
