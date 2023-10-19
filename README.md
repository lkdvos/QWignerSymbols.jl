# QWignerSymbols

Julia package for the implementation of q-deformed Wigner Symbols. Additionally, this
provides an extension to [TensorKit.jl](https://github.com/Jutho/TensorKit.jl) for working
with tensors that have q-deformed SU(2) symmetry.

[![Build Status](https://github.com/lkdvos/QWignerSymbols.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/lkdvos/QWignerSymbols.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/lkdvos/QWignerSymbols.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/lkdvos/QWignerSymbols.jl)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/Q/QWignerSymbols.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/report.html)

Currently, the package provides the following exported functions that define
[q-analogs](https://en.wikipedia.org/wiki/Q-analog):

- `q_number(n::Integer, q::Number)`: $\[n\]_q = \frac{1 - q^n}{1 - q}$
- `q_factorial(n::Integer, q::Number)`: $\[n\]_q! = \prod \[k\]_q$
- `q_binomial(n::Integer, k::Integer, q::Number)`: $\binom{n}{k}_q = \frac{\[n\]_q!}{\[k\]_q! \[n-k\]_q!}$

The following functions are exported for the calculation of q-deformed Wigner Symbols, which
serve a similar function as their
[`WignerSymbols.jl`](https://github.com/Jutho/WignerSymbols.jl) q-less counterparts:

- `q_wigner3j(j1, j2, j3, m1, m2, m3, q)`
- `q_clebschgorda(j1, j2, j3, m1, m2, m3, q)`
- `q_wigner6j(j1, j2, j3, j4, j5, j6, q)`
- `q_racahW(j1, j2, J, j3, J12, J23, q)`

Finally, these can be utilized to construct q-deformed symmetric tensors, by using
`SU2qIrrep{q}` as a drop-in replacement for
[`TensorKit.jl`](https://github.com/Jutho/TensorKit.jl)'s `SU2Irrep`:

```julia
using TensorKit, QWignerSymbols

q = 1.1
j = 1//2
irrep = SU2qIrrep{q}(j)

# Construct a rank-4 tensor with q-deformed SU(2) symmetry
t = TensorMap(randn, Float64, Vect[SU2qIrrep{q}](1//2 => 1)^2 ← Vect[SU2qIrrep{q}](1//2 => 1)^2)
```
