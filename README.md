# BiComplex.jl
Bicomplex numbers in pure Julia.

An implementation of bi-complex numbers following [Luna-Elizarraras et al](https://core.ac.uk/download/pdf/215754454.pdf) [1].

If you want derivatives of functions evaluated at complex locations, I recommend using [DualNumbers.jl](https://github.com/JuliaDiff/DualNumbers.jl)
until [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) enables this feature:

```
julia> using DualNumbers

julia> f(x) = exp(-x^2) # value, f(x)
f (generic function with 1 method)

julia> g(x) = - 2x * exp(-x^2) # derivative, df / dx
g (generic function with 1 method)

julia> z = randn(ComplexF64)
0.46103058785783857 - 0.40369871010374947im

julia> DualNumbers.realpart(f(Dual(1.0 + im, 1))) / f(1.0 + im)
1.0 + 0.0im

julia> DualNumbers.dualpart(f(Dual(1.0 + im, 1))) / g(1.0 + im)
1.0 - 0.0im
```

## References
[1] Luna-Elizarraras, M.E, Shapiro, M., Struppa, D.C., & Vajiac, A. (2012). Bicomplex Numbers and their Elementary Functions. CUBO,
14(2), 61-80. doi: 10.4067/S0719-06462012000200004.
