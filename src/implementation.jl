using AutoHashEquals

@auto_hash_equals struct Bicomplex{T} <: Number
  re::Complex{T}
  im::Complex{T}
end

function Base.convert(::Type{Bicomplex{T}}, x::Real) where T
  return Bicomplex(T(x) + 0im, zero(T) + 0im)
end
function Base.convert(::Type{Bicomplex{T}}, x::Complex) where T
  return Bicomplex(T(x), zero(T))
end
function Base.promote_rule(::Type{Bicomplex{T}}, ::Type{U}) where {T, U<:Real}
  return Bicomplex{promote_type(T, U)}
end

Base.first(b::Bicomplex) = b.re
second(b::Bicomplex) = b.im

const imim  = Bicomplex(0 + 0im, 1 + 0im)

Base.isequal(a::Bicomplex, b::Complex) = (first(a) == first(b)) && (second(a) == second(b))

const e⃗2 = Bicomplex(1 + 0im, 0 + im)
const e⃗ᵀ2 = Bicomplex(1 + 0im, 0 - im)

Base.:/(a::Bicomplex, b::Number) = Bicomplex(first(a) / b, second(a) / b)

function _op(op::T, a::Bicomplex) where {T}
  α, β = idempotentelements(a) 
  return (op(α) * e⃗2 + op(β) * e⃗ᵀ2) / 2
end

const e⃗ = e⃗2 / 2
const e⃗ᵀ = e⃗ᵀ2 / 2

idempotentelements(a::Complex, b::Complex) = (a - im * b, a + im * b)
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

Base.angle(a::Bicomplex) = atan(second(a) / first(a))

function Base.:*(a::Bicomplex, b::Bicomplex)
  return Bicomplex(first(a) * first(b) - second(a) * second(b), first(a) * second(b) + second(a) * first(b))
end
Base.:*(a::Bicomplex, b::Number) = Bicomplex(first(a) * b, second(a) * b)
Base.:*(a::Number, b::Bicomplex) = b * a

#Base.inv(a::Bicomplex) = conj(a) / (first(a)^2 + second(a)^2)
Base.:/(a::Bicomplex, b::Bicomplex) = a * _op(inv, b)#a * inv(b)
Base.:/(a::Number, b::Bicomplex) = a * _op(inv, b)#a * inv(b)

#Base.exp(a::Bicomplex) = exp(first(a)) * Bicomplex(cos(second(a)), sin(second(a)))
#Base.cos(a::Bicomplex) = (exp(imim * a) + exp(-imim * a)) / 2
#Base.sin(a::Bicomplex) = (exp(imim * a) - exp(-imim * a)) / 2 / imim
#Base.tan(a::Bicomplex) = sin(a) / cos(a)
Base.:^(a::Bicomplex, x::AbstractFloat) = _op(z -> z^x, a)
Base.exp(a::Bicomplex) = _op(exp, a)
Base.sqrt(a::Bicomplex) = _op(sqrt, a)
Base.cbrt(a::Bicomplex) = _op(cbrt, a)
Base.cos(a::Bicomplex) = _op(cos, a)
Base.sin(a::Bicomplex) = _op(sin, a)
Base.tan(a::Bicomplex) = _op(tan, a)
Base.log(a::Bicomplex) = _op(log, a)

Base.acos(a::Bicomplex) = _op(acos, a)
Base.asin(a::Bicomplex) = _op(asin, a)
Base.atan(a::Bicomplex) = _op(atan, a)
#Base.acos(a::Bicomplex) = - imim * log(a + sqrt(a^2 - 1))
#Base.asin(a::Bicomplex) = - imim * log(imim * a + sqrt(1 - a^2))
#Base.atan(a::Bicomplex) = asin(a / sqrt(1 + a^2))
#Base.atan(a::Bicomplex) = imim * log((1 + imim * a) / (1 - imim * a))

function derivative(f::T, x::Complex{U}, h=sqrt(eps(U))) where {T,U<:Real}
  return (f(Bicomplex(x, Complex(h, 0)))).im / h
end


