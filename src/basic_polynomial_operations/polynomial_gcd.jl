#############################################################################
#############################################################################
#
# This file implements polynomial GCD 
#                                                                               
#############################################################################
#############################################################################

"""
The extended euclid algorithm for polynomials modulo prime.
"""
function extended_euclid_alg(a::Polynomial, b::Polynomial, prime::Int)
    old_r, r = mod(a, prime), mod(b, prime)
    old_s, s = one(Polynomial), zero(Polynomial)
    old_t, t = zero(Polynomial), one(Polynomial)

    while !iszero(mod(r,prime))
        q = first(divide(old_r, r)(prime))
        old_r, r = r, mod(old_r - q*r, prime)
        old_s, s = s, mod(old_s - q*s, prime)
        old_t, t = t, mod(old_t - q*t, prime)
    end
    g, s, t = old_r, old_s, old_t

    @assert mod(s*a + t*b - g, prime) == 0
    return g, s, t  
end

gcd(a::Polynomial, b::Polynomial, prime::Int) = extended_euclid_alg(a,b,prime) |> first


function extended_euclid_alg(a::PolynomialSparse, b::PolynomialSparse, prime::Int)
    old_j, j = mod(a, prime), mod(b, prime)
    old_k, k = one(PolynomialSparse), zero(PolynomialSparse)
    old_l, l = zero(PolynomialSparse), one(PolynomialSparse)

    while !iszero(mod(j, prime))
        p = first(divide(old_j, j)(prime))
        old_j, j = j, mod(old_j - p*j, prime)
        jre = PolynomialSparse([i for i in j.terms])
        j = jre
        old_k, k = k, mod(old_k - p*k, prime)
        kre = PolynomialSparse([i for i in k.terms])
        k = kre
        old_l, l = l, mod(old_l - p*l, prime)
        lre = PolynomialSparse([i for i in l.terms])
        l = lre
    end
    g, s, t = old_j, old_k, old_l

    @assert mod(s*a + t*b - g, prime) == 0
    return g, s, t  
end

"""
The GCD of two polynomials modulo prime.
"""
gcd(a::PolynomialSparse, b::PolynomialSparse, prime::Int) = extended_euclid_alg(a,b,prime) |> first


function extended_euclid_alg(a::PolynomialSparseBInt, b::PolynomialSparseBInt, prime::BigInt)
    old_j, j = mod(a, prime), mod(b, prime)
    old_k, k = one(PolynomialSparseBInt), zero(PolynomialSparseBInt)
    old_l, l = zero(PolynomialSparseBInt), one(PolynomialSparseBInt)

    while !iszero(mod(j, prime))
        p = first(divide(old_j, j)(prime))
        old_j, j = j, mod(old_j - p*j, prime)
        jre = PolynomialSparseBInt([i for i in j.terms])
        j = jre
        old_k, k = k, mod(old_k - p*k, prime)
        kre = PolynomialSparseBInt([i for i in k.terms])
        k = kre
        old_l, l = l, mod(old_l - p*l, prime)
        lre = PolynomialSparseBInt([i for i in l.terms])
        l = lre
    end
    g, s, t = old_j, old_k, old_l

    @assert mod(s*a + t*b - g, prime) == 0
    return g, s, t  
end

"""
The GCD of two polynomials modulo prime.
"""
gcd(a::PolynomialSparseBInt, b::PolynomialSparseBInt, prime::BigInt) = extended_euclid_alg(a,b,prime) |> first