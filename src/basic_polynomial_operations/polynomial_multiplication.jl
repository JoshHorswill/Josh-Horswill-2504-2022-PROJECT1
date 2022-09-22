#############################################################################
#############################################################################
#
# This file implements polynomial multiplication 
#                                                                               
#############################################################################
#############################################################################

### POLYNOMIALDENSE MULTIPLICATION & POWER OF SPARSE, BINT
"""
Multiply two polynomials.
"""
function *(p1::PolynomialDense, p2::PolynomialDense)::PolynomialDense
    p_out = PolynomialDense()
    for t in p1
        new_summand = (t * p2)
        p_out = p_out + new_summand
    end
    return p_out
end

"""
Power of a polynomial.
"""
function ^(p::PolynomialDense, n::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
    end
    return out
end

function ^(p::PolynomialSparse, n::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
    end
    return out
end

function ^(p::PolynomialSparseBInt, n::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
    end
    return out
end


### MULTIPLICATION OF SPARSE AND SPARSEBINT FOR POLYNOMIALS, TERMS AND INTS
function *(p1::PolynomialSparse, p2::PolynomialSparse)::PolynomialSparse
    p_out = PolynomialSparse()
    for t in p1
        new_summand = (t * p2)
        p_out = p_out + new_summand
    end
    return p_out
end


function *(p1::PolynomialSparseBInt, p2::PolynomialSparseBInt)::PolynomialSparseBInt
    p_out = PolynomialSparseBInt()
    for t in p1
        new_summand = (t * p2)
        p_out = p_out + new_summand
    end
    return p_out
end

*(t::Term, p1::PolynomialModP)::PolynomialModP = iszero(t) ? PolynomialModP() : PolynomialModP(PolynomialSparse(map((pt)->t*pt, p1.p.terms)), p1.prime)
*(p1::PolynomialModP, t::Term)::PolynomialModP = t*p1
*(n::Int, p::PolynomialModP)::PolynomialModP = p*Term(n,0)
*(p::PolynomialModP, n::Int)::PolynomialModP = n*p

### POLYNOMIAL MOD P MULTIPLICATION FOR POLYNOMIALS, TERMS AND INTS 
# ASSERT EQUAL PRIMES & RETURN MOD P
function *(p1::PolynomialModP, p2::PolynomialModP)::PolynomialModP
    @assert p1.prime == p2.prime
    p_out = PolynomialSparse()
    for t in p1.p
        new_summand = (t * p2.p)
        p_out += new_summand
    end
    p_out = mod(p_out, p1.prime)
    return PolynomialModP(p_out, p1.prime)
end


### POWER OF POLYNOMIAL MOD P

function ^(p::PolynomialModP, n::Int)
    poly = p.p    
    polysq = ^(poly, n)
    polysq = mod(polysq, p.prime)
    return PolynomialModP(polysq, p.prime)
end

### REPEATED SQUARED METHOD OF EXPONENTS FOR BINT, SPARSE AND MODP
function repsqr_exp(p::PolynomialSparseBInt, n::Int)
    @assert n < 2^8
    incr = [1, 2, 4, 8, 16, 32, 64, 128]
    bin_vals = digits(n, base=2, pad=8)
    return_val = TermBI(big"1",0)
    for i in 1:8
        if bin_vals[i] != 0
            return_val *= ^(p, incr[i])
        else
        end
    end
    return return_val
end


function repsqr_exp(p::PolynomialSparse, n::Int)
    @assert n < 2^8
    incr = [1, 2, 4, 8, 16, 32, 64, 128]
    bin_vals = digits(n, base=2, pad=8)
    return_val = Term(1,0)

    for i in 1:8
        if bin_vals[i] != 0
            return_val *= ^(p, incr[i])
        else
        end
    end
    return return_val
end

function repsqr_exp(p::PolynomialModP, n::Int)
    @assert n < 2^8
    valu = ^(p.p, n)
    return mod(valu, p.prime)
end



### POWER MOD 1 AND 2 FOR POLYNOMIAL MODP
function pow_mod(p::PolynomialModP, n::Int)
    n < 0 && error("No negative power")
    out = one(p.p)
    for _ in 1:n
        out *= p.p
        out = mod(out, p.prime)
    end
    return out
end

function repsqr_exp(p::PolynomialSparse, n::Int)
    @assert n < 2^8
    incr = [1, 2, 4, 8, 16, 32, 64, 128]
    bin_vals = digits(n, base=2, pad=8)
    return_val = Term(1,0)
    for i in 1:8
        if bin_vals[i] != 0
            return_val *= ^(p, incr[i])
        else
        end
    end
    return return_val
end


function pow_mod2(p::PolynomialModP, n::Int)
    n < 0 && error("No negative power")
    out = one(p.p)
    pal = p.p
    return_val = repsqr_exp(p, n)

    return return_val
end