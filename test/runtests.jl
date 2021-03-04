using BiComplex
using Random, Test

Random.seed!(0)

@testset "Bicomplex" begin
@testset "+ - * /" begin
  a = Bicomplex(randn(ComplexF64), randn(ComplexF64))
  unit = Bicomplex(1+0im, 0+0im)
  @test a - a == Bicomplex(0 + 0im, 0 + 0im)
  @test a + a == 2 * a
  @test isapprox(a / a, unit, rtol=10eps(), atol=10eps())
  @test isapprox((a * a) / a / a, unit, rtol=10eps(), atol=10eps())
  @test jm * jm == -unit
  @test jm^2 == -unit
  @test jm / jm == unit
end

@testset "floating point" begin
  z = Bicomplex(randn(ComplexF64, 2)...)
  @test isapprox(1 / (1 / z), z, rtol=0, atol=eps())
end

@testset "exp sin cos log tan" begin
  a = Bicomplex(randn(ComplexF64), randn(ComplexF64))
  unit = Bicomplex(1+0im, 0+0im)
  @test isapprox(sin(a)^2 + cos(a)^2, unit, rtol=sqrt(eps()), atol=10eps())
  @test isapprox(acos(cos(a)), a, rtol=sqrt(eps()), atol=10eps())
  @test isapprox(cos(acos(a)), a, rtol=sqrt(eps()), atol=10eps())
  @test isapprox(asin(sin(a)), a, rtol=sqrt(eps()), atol=10eps())
  @test isapprox(sin(asin(a)), a, rtol=sqrt(eps()), atol=10eps())
  @test isapprox(tan(atan(a)), a, rtol=sqrt(eps()), atol=10eps())
  @test isapprox(atan(tan(a)), a, rtol=sqrt(eps()), atol=10eps())
  @test isapprox(log(exp(a)), a, rtol=sqrt(eps()), atol=10eps())
  @test isapprox(exp(log(a)), a, rtol=sqrt(eps()), atol=10eps())
end

@testset "derivative" begin
  f(x) = sin(x) * exp(-x^2) / x^3
  g(x) = (cos(x) - 2 * x * sin(x) - 3 * sin(x) / x) / x^3 * exp(-x^2)
  
  for z ∈ (1.0 + im, 1.0 + 0im, 0.0 + im, rand(ComplexF64), -rand(ComplexF64),
           rand(Float64) - im * rand(Float64), -rand(Float64) + im * rand(Float64))
    @test isapprox(BiComplex.derivative(f, z), g(z), rtol=sqrt(eps()), atol=10eps())
    fop(x) = BiComplex._op(f, x)
    @test isapprox(BiComplex.derivative(fop, z), g(z), rtol=sqrt(eps()), atol=10eps())
    @test isapprox(BiComplex.derivative(exp, z), exp(z), rtol=sqrt(eps()), atol=10eps())
    @test isapprox(BiComplex.derivative(cos, z),-sin(z), rtol=2sqrt(eps()),
                   atol=10eps()) # annoyingly needs rtol of 2 sqrt eps
    @test isapprox(BiComplex.derivative(sin, z), cos(z), rtol=sqrt(eps()), atol=10eps())
    @test isapprox(BiComplex.derivative(tan, z), 1/cos(z)^2, rtol=sqrt(eps()), atol=10eps())
    if (1-z^2) != 0
      @test isapprox(BiComplex.derivative(acos, z),-1/sqrt(1-z^2), rtol=sqrt(eps()), atol=10eps())
      @test isapprox(BiComplex.derivative(asin, z), 1/sqrt(1-z^2), rtol=sqrt(eps()), atol=10eps())
    end
    if (1+z^2) != 0
      @test isapprox(BiComplex.derivative(atan, z), 1/(1+z^2), rtol=sqrt(eps()), atol=10eps())
    end
  end
end

@testset "Promotion conversion, etc" begin
  @test Bicomplex(0.0 + 0im, 1 + 0im) == jm
  @test jm^2 == -1
end

end

