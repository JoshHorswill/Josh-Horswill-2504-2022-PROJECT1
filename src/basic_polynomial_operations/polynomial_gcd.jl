#############################################################################
#############################################################################
#
# This file implements polynomial GCD 
#                                                                               
#############################################################################
#############################################################################

### EUCLIDEAN ALGORITHM AND GCD DEFINITIONS FOR DENSE, SPARSE & BINT
function extended_euclid_alg(a::PolynomialDense, b::PolynomialDense, prime::Int)
    old_r, r = mod(a, prime), mod(b, prime)
    old_s, s = one(PolynomialDense), zero(PolynomialDense)
    old_t, t = zero(PolynomialDense), one(PolynomialDense)

    while !iszero(mod(r,prime))
        q = first(divide(old_r, r)(prime))
        old_r, r = r, mod(old_r - q*r, prime)
        old_s, s = s, mod(old_s - q*s, prime)
        old_t, t = t, mod(old_t - q*t, prime)
    end
    g, s, t = old_r, old_s, old_t

    plswork = mod(s*a + t*b - g, prime)
    @assert degree(plswork) == 0
    return g, s, t  
end

gcd(a::PolynomialDense, b::PolynomialDense, prime::Int) = extended_euclid_alg(a,b,prime) |> first


function extended_euclid_alg(a::PolynomialSparse, b::PolynomialSparse, prime::Int)
    old_j, j = mod(a, prime), mod(b, prime)
    old_k, k = one(PolynomialSparse), zero(PolynomialSparse)
    old_l, l = zero(PolynomialSparse), one(PolynomialSparse)

    while !iszero(mod(j, prime))
        p = first(divide(PolynomialSparse([i for i in old_j.terms]), j)(prime))
        old_j, j = j, mod(old_j - p*j, prime)
        jre = PolynomialSparse([i for i in j.terms]);
        j = jre;
        old_k, k = k, mod(old_k - p*k, prime)
        kre = PolynomialSparse([i for i in k.terms]);
        k = kre;
        old_l, l = l, mod(old_l - p*l, prime)
        lre = PolynomialSparse([i for i in l.terms]);
        l = lre;        
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

    while !iszero2(mod(j, prime))
        p = first(divide(PolynomialSparseBInt([i for i in old_j.terms]), j)(prime))
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

    xx = mod((s*a + t*b) - g, prime)
    degrees(xx)[1] == 0 && coeffs(xx)[1] == 0
    
    return g, s, t  
end


"""
The GCD of two polynomials modulo prime.
"""
gcd(a::PolynomialSparseBInt, b::PolynomialSparseBInt, prime::BigInt) = extended_euclid_alg(a,b,prime) |> first