#############################################################################
#############################################################################
#
# This file defines the TermBI type with several operations 
#                                                                               
#############################################################################
#############################################################################

##############################
# Term type and construction #
##############################

"""
A term for BigInt
"""
struct TermBI  #structs are immutable by default
    coeff::BigInt 
    degree::Int
    function TermBI(coeff::BigInt, degree::Int)
        degree < 0 && error("Degree must be non-negative")
        coeff != 0 ? new(coeff,degree) : new(coeff,0)
    end
end

"""
Creates the zero term.
"""
zero(::Type{TermBI})::TermBI = TermBI(big"0",0)
"""
Creates the unit term.
"""
one(::Type{TermBI})::TermBI = TermBI(big"1",0)

########################################################################################
# Display #
########################################################################################
"""
Show a term.
"""
show(io::IO, t::TermBI) = print(io, "$(t.coeff)⋅x^$(t.degree)") #\cdot + [TAB]


########################################################################################
# Queries about a term #
########################################################################################

"""
Check if a term is 0.
"""
iszero(t::TermBI)::Bool = iszero(t.coeff)
iszero2(t::TermBI)::Bool = iszero(t.degree)

"""
Compare two terms.
"""
isless(t1::TermBI,t2::TermBI)::Bool =  t1.degree == t2.degree ? (t1.coeff < t2.coeff) : (t1.degree < t2.degree)

"""
Evaluate a term at a point x.
"""
evaluate(t::TermBI, x::T) where T <: Number = t.coeff * x^t.degree

########################################################################################
# Operations with a term #
########################################################################################

"""
Add two terms of the same degree.
"""
function +(t1::TermBI,t2::TermBI)::TermBI
    @assert t1.degree == t2.degree
    TermBI(t1.coeff + t2.coeff, t1.degree)
end

"""
Negate a term.
"""
-(t::TermBI,) = TermBI(-t.coeff,t.degree)

"""
Subtract two terms with the same degree.
"""
-(t1::TermBI, t2::TermBI)::TermBI = t1 + (-t2)
"""
Multiply two terms.
"""
*(t1::TermBI, t2::TermBI)::TermBI = TermBI(t1.coeff * t2.coeff, t1.degree + t2.degree)
#*(t1::TermBI, t2::Int64)::TermBI = TermBI(t1.coeff * t2, t1.degree)
"""
Compute the symmetric mod of a term with an integer.
"""
mod(t::TermBI, p::BigInt) = TermBI(mod(t.coeff,p), t.degree)

"""
Compute the derivative of a term.
"""
derivative(t::TermBI) = TermBI(t.coeff*t.degree,max(t.degree-1,0))

"""
Divide two terms. Returns a function of an integer.
"""
function ÷(t1::TermBI,t2::TermBI) #\div + [TAB]
    @assert t1.degree ≥ t2.degree
    f(p::BigInt)::TermBI = TermBI(mod((t1.coeff * int_inverse_mod(t2.coeff, p)), p), t1.degree - t2.degree)
end

"""
Integer divide a term by an integer.
"""
÷(t::TermBI, n::BigInt) = t ÷ TermBI(n,0)

q = 2