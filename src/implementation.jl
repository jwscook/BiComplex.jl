using AutoHashEquals

@auto_hash_equals struct Bicomplex{T} <: Number
  re::Complex{T}
  im::Complex{T}
end
Bicomplex(a::Complex{T}, b::Complex{U}) where {T,U} = Bicomplex(promote(a, b)...)
function Bicomplex{T}(a::Bicomplex{U}) where {T,U}
  V = complex(promote_type(T, U))
  return Bicomplex(V(first(a)), V(second(a)))
end

function Base.convert(::Type{Bicomplex{T}}, x::Bicomplex{U}) where {T, U}
  V = complex(promote_type(T, U))
  return Bicomplex(V(x.re), V(x.im))
end

function Base.convert(::Type{Bicomplex{T}}, x::U) where {T, U<:Number}
  V = complex(promote_type(T, U))
  return Bicomplex(V(x), V(0))
end

const jm = Bicomplex(Complex(false, false), Complex(true, false))

Base.first(b::Bicomplex) = b.re
second(b::Bicomplex) = b.im

function Base.promote_rule(::Type{Bicomplex{T}}, ::Type{Bicomplex{U}}
    ) where {T<:Real, U<:Real}
  return Bicomplex{promote_type(T, U)}
end
function Base.promote_rule(::Type{Bicomplex{T}}, ::Type{U}
    ) where {T<:Real, U<:Number}
  return Bicomplex{promote_type(T, real(U))}
end

Base.isfinite(a::Bicomplex) = isfinite(first(a)) && isfinite(second(a))
Base.isnan(a::Bicomplex) = isnan(first(a)) || isnan(second(a))

const e⃗2 = Bicomplex(1 + 0im, 0 + im)
const e⃗ᵀ2 = Bicomplex(1 + 0im, 0 - im)

Base.:/(a::Bicomplex, b::Number) = Bicomplex(first(a) / b, second(a) / b)

function _op(op::T, a::Bicomplex) where {T}
  α, β = idempotentelements(a) 
  return (op(α) * e⃗2 + op(β) * e⃗ᵀ2) / 2
end

const e⃗ = e⃗2 / 2
const e⃗ᵀ = e⃗ᵀ2 / 2

idempotentelements(a::Complex, b::Complex) = (imb = im * b; (a - imb, a + imb))
idempotentelements(a::Bicomplex) = idempotentelements(first(a), second(a))

Base.conj(a::Bicomplex) = Bicomplex(first(a), -second(a))
Base.abs2(a::Bicomplex) = abs2(first(a)) + abs2(second(a))
Base.abs(a::Bicomplex) = sqrt(abs2(first(a)) + abs2(second(a)))
Base.rand(::Type{Bicomplex{T}}) where {T} = Bicomplex(rand(T), rand(T))
Base.randn(::Type{Bicomplex{T}}) where {T} = Bicomplex(randn(T), randn(T))

Base.:+(a::Bicomplex, b::Bicomplex) = Bicomplex(first(a) + first(b), second(a) + second(b))
Base.:-(a::Bicomplex, b::Bicomplex) = Bicomplex(first(a) - first(b), second(a) - second(b))
Base.:+(a::Bicomplex, b::Number) = Bicomplex(first(a) + b, second(a))
Base.:-(a::Bicomplex, b::Number) = Bicomplex(first(a) - b, second(a))
Base.:-(a::Number, b::Bicomplex) = Bicomplex(a - first(b), -second(b))
Base.:-(a::Bicomplex) = Bicomplex(-first(a), -second(a))

Base.:angle(a::Bicomplex) = atan(second(a) / first(a))

function Base.:*(a::Bicomplex, b::Bicomplex)
  return Bicomplex(first(a) * first(b) - second(a) * second(b), first(a) * second(b) + second(a) * first(b))
end
Base.:*(a::Bicomplex, b::Number) = Bicomplex(first(a) * b, second(a) * b)
Base.:*(a::Number, b::Bicomplex) = b * a

Base.:/(a::Bicomplex, b::Bicomplex) = a * _op(inv, b)
Base.:/(a::Number, b::Bicomplex) = a * _op(inv, b)

Base.:^(a::Bicomplex, x::AbstractFloat) = _op(z -> z^x, a)
for op ∈ (:sqrt, :cbrt, :exp, :log, :cos, :sin, :tan, :acos, :asin, :atan)
  @eval Base.$op(x::Bicomplex) = _op($op, x)
end

function derivative(f::T, x::U, h=sqrt(eps(real(U)))) where {T,U<:Number}
  return U(valuederivative(f, x, h)[2])
end

function valuederivative(f::T, x::U, h=sqrt(eps(real(U)))) where {T,U<:Number}
  fx = f(Bicomplex(Complex(x), Complex(h)))
  value = U(first(fx))
  deriv = U(second(fx)) / h
  return (value, deriv)
end

function valuederivatives(f::T, x::U, h=cbrt(eps(U))) where {T,U<:Real}
  fx = f(Bicomplex(Complex(x, h), Complex(h, 0)))
  value = real(first(fx))
  firstderiv = imag(first(fx)) / h
  secondderiv = imag(second(fx)) / h^2
  return (value, firstderiv, secondderiv)
end

function secondderivative(f::T, x::U, h=cbrt(eps(U))) where {T,U<:Real}
  return valuederivatives(f, x, h)[3]
end

