#############################################################################
#############################################################################
#
# This file implements polynomial addition 
#                                                                               
#############################################################################
#############################################################################

### ADDITION OF POLYNOMIAL TYPES DENSE, SPARSE & BINT, W.R.T TO POLYNOMIALS & TERMS
function +(p::PolynomialDense, t::Term)
    p = deepcopy(p)
    if t.degree > degree(p)
        push!(p, t)
    else
        if !iszero(p.terms[t.degree + 1]) #+1 is due to indexing
            p.terms[t.degree + 1] += t
        else
            p.terms[t.degree + 1] = t
        end
    end

    return trim!(p)
end
+(t::Term, p::PolynomialDense) = p + t

function +(p::PolynomialSparse, t::Term)
    p = deepcopy(p)
    cnt = 0

    for terms in p
        cnt += 1
        if terms.degree == t.degree
            p.terms[cnt] += t
        else
            p.terms[cnt] = terms
        end
    end

    if t.degree ∉ degrees(p)
        push!(p, t)
    end
    return p
end

function +(p::PolynomialSparseBInt, t::TermBI)
    p = deepcopy(p)
    cnt = 0

    for terms in p
        cnt += 1
        if terms.degree == t.degree
            p.terms[cnt] += t
        else
            p.terms[cnt] = terms
        end
    end

    if t.degree ∉ degrees(p)
        push!(p, t)
    end
    return p
end

+(t::Term, p::PolynomialSparse) = p + t
+(t::TermBI, p::PolynomialSparseBInt) = p + t

"""
Add two polynomials.
"""
function +(p1::PolynomialDense, p2::PolynomialDense)::PolynomialDense
    p = deepcopy(p1)
    for t in p2
        p += t
    end
    return p
end

function +(p1::PolynomialSparse, p2::PolynomialSparse)::PolynomialSparse
    p = deepcopy(p1)
    for t in p2
        p += t
    end
    return p
end

function +(p1::PolynomialSparseBInt, p2::PolynomialSparseBInt)::PolynomialSparseBInt
    p = deepcopy(p1)
    for t in p2
        p += t
    end
    return p
end


### POLYNOMIAL MOD P IMPLEMENTATION FOR TERMS, INTS AND POLYNOMIALS
function +(f::PolynomialModP, t::Term)
    p = deepcopy(f.p)

    cnt = 0
    for terms in p
        cnt += 1
        if terms.degree == t.degree
            p.terms[cnt] += t
        else 
            p.terms[cnt] = terms
        end
    end

    if t.degree ∉ degrees(f)
        push!(p, t)
    end
    return PolynomialModP(p, f.prime)
end

+(t::Term, p::PolynomialModP) = p + t

function +(f1::PolynomialModP, f2::PolynomialModP)::PolynomialModP
    @assert f1.prime == f2.prime
    p = deepcopy(f1.p)
    for t in f2.p
        p += t
    end
    p_out = mod(p, f1.prime)
    return PolynomialModP(p_out, f1.prime)
end

+(p::PolynomialModP, n::Int) = p + Term(n,0)
+(n::Int, p::PolynomialModP) = p + Term(n,0)
-(p::PolynomialModP, t::Term) = p + (-t)



### GENERAL NOTATION
+(p::PolynomialDense, n::Int) = p + Term(n,0)
+(n::Int, p::PolynomialDense) = p + Term(n,0)

+(p::PolynomialSparse, n::Int) = p + Term(n,0)
+(n::Int, p::PolynomialSparse) = p + Term(n,0)

+(p::PolynomialSparseBInt, n::BigInt) = p + TermBI(n,0)
+(n::BigInt, p::PolynomialSparseBInt) = p + TermBI(n,0)

-(p::PolynomialSparse, t::Term) = p + (-t)
-(p::PolynomialSparseBInt, t::TermBI) = p + (-t)