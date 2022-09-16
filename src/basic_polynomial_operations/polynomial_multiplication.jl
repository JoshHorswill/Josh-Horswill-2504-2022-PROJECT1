#############################################################################
#############################################################################
#
# This file implements polynomial multiplication 
#                                                                               
#############################################################################
#############################################################################

"""
Multiply two polynomials.
"""
function *(p1::Polynomial, p2::Polynomial)::Polynomial
    p_out = Polynomial()
    for t in p1
        new_summand = (t * p2)
        p_out = p_out + new_summand
    end
    return p_out
end

"""
Power of a polynomial.
"""
function ^(p::Polynomial, n::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
    end
    return out
end

function ^(p::PolynomialSparseBIn, n::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
    end
    return out
end

function ^(p::PolynomialSparse, n::Int)::PolynomialSparse

    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
    end

    println(length(out))
    terms = [zero(Term) for _ in 0:100]
    cnter = 0

    for t in out
        cnter += 1
        #println(t)
        if t.coeff !=0 && t.degree != 0 
            cnter += 1
            terms[cnter] = Term(t.coeff, t.degree-1) #+1 accounts for 1-indexing
        end
        #println(terms)
    end
    
    cnter_2 = 0
    return_val = PolynomialSparse(Term(0,0))
    for i in terms
        cnter_2 += 1
        if i.degree !== 0 && i.coeff !== 0
            return_val += i
        end
    end
    return return_val
end

function ^(p::PolynomialSparseBIn, n::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
    end

    len = length(out)
    terms = [zero(TermBI) for _ in 0:100]
    println(typeof(terms))
    cnter = 0

    for t in out
        cnter += 1
        if t.coeff !=0 && t.degree != 0 
            cnter += 1
            terms[cnter] = TermBI(t.coeff, t.degree-1) #+1 accounts for 1-indexing
        end
    end

    cnter_2 = 0
    return_val = PolynomialSparseBIn(TermBI(big"0",0))
    for i in terms
        cnter_2 += 1
        if i.degree !== 0 && i.coeff !== 0
            return_val += i
        end
    end

    return return_val
end

function ^(p::PolynomialSparseBInt, n::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
    end

    len = length(out)
    terms = [zero(TermBI) for _ in 0:100]
    println(typeof(terms))
    cnter = 0

    for t in out
        cnter += 1
        if t.coeff !=0 && t.degree != 0 
            cnter += 1
            terms[cnter] = TermBI(t.coeff, t.degree-1) #+1 accounts for 1-indexing
        end
    end

    cnter_2 = 0
    return_val = PolynomialSparseBInt(TermBI(big"0",0))
    for i in terms
        cnter_2 += 1
        if i.degree !== 0 && i.coeff !== 0
            return_val += i
        end
    end

    return return_val
end

function *(p1::PolynomialSparse, p2::PolynomialSparse)::PolynomialSparse
    p_out = PolynomialSparse()
    for t in p1
        new_summand = (t * p2)
        p_out = p_out + new_summand
    end
    return p_out
end

function *(p1::PolynomialSparseBIn, p2::PolynomialSparseBIn)::PolynomialSparseBIn
    p_out = PolynomialSparseBIn()
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


"""
Power of a polynomial.
"""

