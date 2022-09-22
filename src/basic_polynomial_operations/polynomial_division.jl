#############################################################################
#############################################################################
#
# This file implements polynomial division 
#                                                                               
#############################################################################
#############################################################################

### DEFINITION OF DIVISION OF POLYNOMIAL DENSE AND NUM, DEM OPERATORS
function divide(num::PolynomialDense, den::PolynomialDense)
    function division_function(p::Int)
        f, g = mod(num,p), mod(den,p)
        degree(f) < degree(num) && return nothing 
        iszero(g) && throw(DivideError())
        q = PolynomialDense()
        prev_degree = degree(f)
        while degree(f) ≥ degree(g) 
            h = PolynomialDense( (leading(f) ÷ leading(g))(p) )  #syzergy 
            f = mod((f - h*g), p)
            q = mod((q + h), p)  
            prev_degree == degree(f) && break
            prev_degree = degree(f)
        end
        @assert iszero( mod((num  - (q*g + f)),p))
        return q, f
    end
    return division_function
end

"""
The quotient from polynomial division. Returns a function of an integer.
"""
÷(num::PolynomialDense, den::PolynomialDense)  = (p::Int) -> first(divide(num,den)(p))
"""
The remainder from polynomial division. Returns a function of an integer.
"""
rem(num::PolynomialDense, den::PolynomialDense)  = (p::Int) -> last(divide(num,den)(p))



### IMPLEMENTATION OF DIVISION FOR POLYNOMIAL SPARSE AND SUBSEQUENT OPERATORS
# Due to ongoing issue, repacking the values was implemented, increases computation slightly but completely effective.
function divide(num::PolynomialSparse, den::PolynomialSparse)
    function division_function(p::Int)
        f, g = mod(num,p), mod(den,p)
        degree(f) < degree(num) && return nothing 
        iszero(g) && throw(DivideError())
        q = Term(0,0)
        prev_degree = degree(f)
        while degree(f) ≥ degree(g) 
            h = PolynomialSparse( (leading(f) ÷ leading(g))(p) )  #syzergy 
            f = mod((f - h*g), p)
            fre = PolynomialSparse([i for i in f.terms])
            f = fre
            q = mod((q + h), p)  
            qre = PolynomialSparse([i for i in q.terms])
            q = qre
            prev_degree == degree(f) && break
            prev_degree = degree(f)
        end
        @assert iszero( mod((num  - (q*g + f)),p))
        return q, f
    end
    return division_function
end
"""
The quotient from polynomial division. Returns a function of an integer.
"""
÷(num::PolynomialSparse, den::PolynomialSparse)  = (p::Int) -> first(divide(num,den)(p))
"""
The remainder from polynomial division. Returns a function of an integer.
"""
rem(num::PolynomialSparse, den::PolynomialSparse)  = (p::Int) -> last(divide(num,den)(p))



### POLYNOMIAL SPARSE BINT IMPLEMENTATION
# Repacking used again, follows on from the sparse interpretation just with TermBI

function divide(num::PolynomialSparseBInt, den::PolynomialSparseBInt)
    function division_function(p::BigInt)
        f, g = mod(num,p), mod(den,p)
        degree(f) < degree(num) && return nothing 
        iszero(g) && throw(DivideError())
        q = TermBI(big"0",0)
        prev_degree = degree(f)
        while degree(f) ≥ degree(g) 
            h = PolynomialSparseBInt( (leading(f) ÷ leading(g))(p) )  #syzergy 
            f = mod((f - h*g), p)
            fre = PolynomialSparseBInt([i for i in f.terms])
            f = fre
            q = mod((q + h), p)  
            qre = PolynomialSparseBInt([i for i in q.terms])
            q = qre
            prev_degree == degree(f) && break
            prev_degree = degree(f)
        end
        @assert iszero2( mod((num  - (q*g + f)),p))
        return q, f
    end
    return division_function
end
 

"""
The quotient from polynomial division. Returns a function of an integer.
"""
÷(num::PolynomialSparseBInt, den::PolynomialSparseBInt)  = (p::Int) -> first(divide(num,den)(p))
"""
The remainder from polynomial division. Returns a function of an integer.
"""
rem(num::PolynomialSparseBInt, den::PolynomialSparseBInt)  = (p::Int) -> last(divide(num,den)(p))


### POLYNOMIAL MOD P IMPLEMENTATION OF DIVISION
# ASSERT PRIMES ARE THE SAME, RETURN IN MOD P
function divide(num::PolynomialModP, den::PolynomialModP)
    @assert num.prime == den.prime

    num_val = num.p
    den_val = den.p
    p = num.prime

    function division_function(p::Int)
        f, g = mod(num_val,p), mod(den_val,p)
        degree(f) < degree(num_val) && return nothing 
        iszero(g) && throw(DivideError())
        q = Term(0,0)
        prev_degree = degree(f)
        while degree(f) ≥ degree(g) 
            h = PolynomialSparse( (leading(f) ÷ leading(g))(p) )  #syzergy 
            f = mod((f - h*g), p)
            fre = PolynomialSparse([i for i in f.terms])
            f = fre
            q = mod((q + h), p)  
            qre = PolynomialSparse([i for i in q.terms])
            q = qre
            prev_degree == degree(f) && break
            prev_degree = degree(f)
        end
        @assert iszero( mod((num_val  - (q*g + f)),p))
        return q, f
    end

    return division_function(p)
end


function ^(p::PolynomialModP, n::Int)
    poly = p.p    
    polysq = ^(poly, n)
    return PolynomialModP(polysq, p.prime)
end