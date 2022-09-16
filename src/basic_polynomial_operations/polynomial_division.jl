#############################################################################
#############################################################################
#
# This file implements polynomial division 
#                                                                               
#############################################################################
#############################################################################

"""  Modular algorithm.
f divide by g

f = q*g + r

p is a prime
"""
function divide(num::Polynomial, den::Polynomial)
    function division_function(p::Int)
        f, g = mod(num,p), mod(den,p)
        degree(f) < degree(num) && return nothing 
        iszero(g) && throw(DivideError())
        q = Polynomial()
        prev_degree = degree(f)
        while degree(f) ≥ degree(g) 
            h = Polynomial( (leading(f) ÷ leading(g))(p) )  #syzergy 
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
÷(num::Polynomial, den::Polynomial)  = (p::Int) -> first(divide(num,den)(p))
"""
The remainder from polynomial division. Returns a function of an integer.
"""
rem(num::Polynomial, den::Polynomial)  = (p::Int) -> last(divide(num,den)(p))



# POLYNOMIAL SPARSE
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


# POLYNOMIAL SPARSE BI
function divide(num::PolynomialSparseBIn, den::PolynomialSparseBIn)
    function division_function(p::Int)
        f, g = mod(num,p), mod(den,p)
        degree(f) < degree(num) && return nothing 
        #iszero(g) && throw(DivideError())
        q = TermBI(big"0",0)
        prev_degree = degree(f)
        while degree(f) ≥ degree(g) 
            h = PolynomialSparseBIn( (leading(f) ÷ leading(g))(p) )  #syzergy 
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
÷(num::PolynomialSparseBIn, den::PolynomialSparseBIn)  = (p::Int) -> first(divide(num,den)(p))

"""
The remainder from polynomial division. Returns a function of an integer.
"""
rem(num::PolynomialSparseBIn, den::PolynomialSparseBIn)  = (p::Int) -> last(divide(num,den)(p))

# POLYNOMIAL SPARSE
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