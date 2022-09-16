#############################################################################
#############################################################################
#
# This file implements polynomial addition 
#                                                                               
#############################################################################
#############################################################################

"""
Add a polynomial and a term.
"""
function +(p::Polynomial, t::Term)
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
+(t::Term, p::Polynomial) = p + t

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

function +(p::PolynomialSparseBIn, t::TermBI)
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
+(t::TermBI, p::PolynomialSparseBIn) = p + t
+(t::TermBI, p::PolynomialSparseBInt) = p + t
#+(p::PolynomialSparseBIn, t::TermBI) = p + t

"""
Add two polynomials.
"""
function +(p1::Polynomial, p2::Polynomial)::Polynomial
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

function +(p1::PolynomialSparseBIn, p2::PolynomialSparseBIn)::PolynomialSparseBIn
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

"""
Add a polynomial and an integer.
"""
+(p::Polynomial, n::Int) = p + Term(n,0)
+(n::Int, p::Polynomial) = p + Term(n,0)

+(p::PolynomialSparse, n::Int) = p + Term(n,0)
+(n::Int, p::PolynomialSparse) = p + Term(n,0)

+(p::PolynomialSparseBIn, n::Int) = p + TermBI(n,0)
+(n::Int, p::PolynomialSparseBIn) = p + TermBI(n,0)

+(p::PolynomialSparseBInt, n::BigInt) = p + TermBI(n,0)
+(n::BigInt, p::PolynomialSparseBInt) = p + TermBI(n,0)

-(p::PolynomialSparse, t::Term) = p + (-t)
-(p::PolynomialSparseBIn, t::TermBI) = p + (-t)
-(p::PolynomialSparseBInt, t::TermBI) = p + (-t)