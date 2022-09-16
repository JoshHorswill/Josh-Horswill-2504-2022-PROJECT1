#############################################################################
#############################################################################
#
# This is the main project file for polynomial factorization
#                                                                               
#############################################################################
#############################################################################

using Distributions, StatsBase, Random

import Base: %
import Base: push!, pop!, iszero, show, isless, map, map!, iterate, length, last
import Base: +, -, *, mod, %, ÷, ==, ^, rand, rem, zero, one

include("src/general_alg.jl")
include("src/TermTypes/term.jl")
include("src/TermTypes/termBI.jl")
include("src/PolynomialTypes/polynomial.jl")
include("src/PolynomialTypes/PolynomialSparse.jl")
include("src/PolynomialTypes/PolynomialSparseBI.jl")
include("src/PolynomialTypes/PolynomialModP.jl")
include("src/PolynomialTypes/BigIntPS.jl")
    include("src/basic_polynomial_operations/polynomial_addition.jl")
    include("src/basic_polynomial_operations/polynomial_multiplication.jl")
    include("src/basic_polynomial_operations/polynomial_division.jl")
    include("src/basic_polynomial_operations/polynomial_gcd.jl")
include("src/polynomial_factorization/factor.jl")

println("done")
