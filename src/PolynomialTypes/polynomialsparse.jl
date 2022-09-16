"""
A Polynomial type - designed to be for polynomials with integer coefficients.
"""

struct PolynomialSparse
    
    terms::Vector{Term}   
    
    #Inner constructor of 0 polynomial
    #PolynomialSparse() = new([Term(0,0)])
    PolynomialSparse() = new([zero(Term)])
    
    #Inner constructor of polynomial based on arbitrary list of terms
    function PolynomialSparse(vt::Vector{Term})

        #Filter the vector so that there is not more than a single zero term
        vt = filter((t)->!iszero(t), vt)
        if isempty(vt)
            vt = [zero(Term)]
        end

        # Take count of elements in vector of terms
        cnt = 0
        for t in vt
            cnt += 1
        end
        terms = [zero(Term) for i in 0:cnt-1] #Array of zeroes to the size of the count

        # For each term, if not a null item, commit to terms array and icrement index
        cnter = 0
        inner_ct = 0
        for t in vt
            inner_ct += 1
            if t.coeff !=0 && t.degree != 0
                cnter += 1
                terms[cnter] = t #+1 accounts for 1-indexing
            elseif t.coeff != 0 && t.degree == 0
                cnter += 1
                terms[cnter] = t
            end
        end

        # I am disgusted in myself for having to use this, however here I am.
        ## bubble_poly_sort - sorts sparse polynomials terms in order of reverse degree
        function bubble_poly_sort!(terms)
            n = length(terms)
            deg = [t.degree for t in terms]
            for i in 1:n-1
                for j in 1:n-i
                    if deg[j] < deg[j+1]
                        deg[j], deg[j+1] = deg[j+1], deg[j]
                        terms[j], terms[j+1] = terms[j+1], terms[j]
                    end
                end
            end
            return terms
        end

        bubble_poly_sort!(terms)

        return new(terms)
    end
end

"""
This function maintains the invariant of the Polynomial type so that there are no zero terms beyond the highest
non-zero term.
"""
function trim!(p::PolynomialSparse)::PolynomialSparse
    i = length(p.terms)
    while i > 1
        if iszero(p.terms[i])
            pop!(p.terms)
        else
            break
        end
        i -= 1
    end
    return p
end

"""
Construct a polynomial with a single term.
"""
PolynomialSparse(t::Term) = PolynomialSparse([t])

"""
Construct a polynomial of the form x^p-x.
"""
cyclotonic_polynomial(p::Int) = PolynomialSparse([Term(1,p), Term(-1,0)])


"""
Construct a polynomial of the form x-n.
"""
linear_monic_polynomial(n::Int) = PolynomialSparse([Term(1,1), Term(-n,0)])

"""
Construct a polynomial of the form x.
"""
x_poly_sp() = PolynomialSparse(Term(1,1))

"""
Creates the zero polynomial.
"""
#zero(::Type{PolynomialSparse})::PolynomialSparse = PolynomialSparse()
zero(::Type{PolynomialSparse})::PolynomialSparse = PolynomialSparse(Term(0,0))
"""
Creates the unit polynomial.
"""
#one(::Type{PolynomialSparse})::PolynomialSparse = PolynomialSparse(one(Term))
#one(p::PolynomialSparse) = one(typeof(p))
one(::Type{PolynomialSparse})::PolynomialSparse = PolynomialSparse(Term(1,0))
one(p::PolynomialSparse)::PolynomialSparse = PolynomialSparse(Term(1,0))
"""
Generates a random polynomial.
"""
function rand(::Type{PolynomialSparse} ; 
                degree::Int = -1, 
                terms::Int = -1, 
                max_coeff::Int = 100, 
                mean_degree::Float64 = 5.0,
                prob_term::Float64  = 0.7,
                monic = false,
                condition = (p)->true)
        
    while true 
        _degree = degree == -1 ? rand(Poisson(mean_degree)) : degree
        _terms = terms == -1 ? rand(Binomial(_degree,prob_term)) : terms
        degrees = vcat(sort(sample(0:_degree-1,_terms,replace = false)),_degree)
        coeffs = rand(1:max_coeff,_terms+1)
        monic && (coeffs[end] = 1)
        p = PolynomialSparse( [Term(coeffs[i],degrees[i]) for i in 1:length(degrees)] )
        condition(p) && return p
    end
end

###########
# Display #
###########

"""
Show a polynomial.
"""
function show(io::IO, p::PolynomialSparse) 
    if iszero(p)
        print(io,"0")
    else
        n = length(p.terms)
        for (i,t) in enumerate(p.terms)
            if !iszero(t)
                if t.degree == 0 && t.coeff < 0
                    print(io, i != 1 ? " " : "", t.coeff, " ")
                elseif t.degree == 0 && t.coeff > 0
                    print(io, i != 1 ? " + " : "", t.coeff, "⋅x", "")
                elseif t.degree > 0 && t.coeff > 0
                    print(io, i != 1 ? " + " : "", t)
                elseif t.degree > 0 && t.coeff < 0
                    print(io, i != 1 ? " " : "", t)
                else
                    #print("io, t, i != n ? " + " : """)
                end
            end
        end
    end
end

##############################################
# Iteration over the terms of the polynomial #
##############################################

"""
Allows to do iteration over the non-zero terms of the polynomial. This implements the iteration interface.
"""
iterate(p::PolynomialSparse, state=1) = iterate(p.terms, state)

##############################
# Queries about a polynomial #
##############################

"""
The number of terms of the polynomial.
"""
length(p::PolynomialSparse) = length(p.terms) 

"""
The leading term of the polynomial.
"""
#leading(p::PolynomialSparse)::Term = isempty(p.terms) ? zero(Term) : first(p.terms) 
function leading(p::PolynomialSparse)::Term
    lead = -Inf
    cnt = 1
    co, de = 0, 0
    for terms in p
        #cnt += 1
        if terms.degree > lead
            lead = terms.degree
            co = terms.coeff
            de = terms.degree
        end
    end
    return Term(co, de)
end

"""
Returns the coefficients of the polynomial.
"""
coeffs(p::PolynomialSparse)::Vector{Int} = [t.coeff for t in p]

"""
The degree of the polynomial.
"""
#degree(p::PolynomialSparse)::Int = leading(p).degree 
degrees(p::PolynomialSparse)::Vector{Int} = [t.degree for t in p]

function degree(p::PolynomialSparse)::Int
    max_deg = -Inf
    for terms in p
        if terms.degree > max_deg
            max_deg = terms.degree
        end
    end
    return max_deg
end
"""
The content of the polynomial is the GCD of its coefficients.
"""
content(p::PolynomialSparse)::Int = euclid_alg(coeffs(p))

"""
Evaluate the polynomial at a point `x`.
"""
evaluate(f::Polynomial, x::T) where T <: Number = sum(evaluate(t,x) for t in f)

################################
# Pushing and popping of terms #
################################

"""
Push a new term into the polynomial.
"""
#Note that ideally this would throw and error if pushing another term of degree that is already in the polynomial
function push!(p::PolynomialSparse, t::Term) 
    #if t.degree <= degree(p)
    #    p.terms[1] = t
    #else
    #    append!(p.terms, zeros(Term, t.degree - degree(p)-1))
    push!(p.terms, t)
    #end
    #return p        
end

"""
Pop the leading term out of the polynomial. When polynomial is 0, keep popping out 0.
"""
function pop!(p::PolynomialSparse)::Term 
    popped_term = pop!(p.terms) #last element popped is leading coefficient

    while !isempty(p.terms) && iszero(last(p.terms))
        pop!(p.terms)
    end

    if isempty(p.terms)
        push!(p.terms, zero(Term))
    end

    return popped_term
end

"""
Check if the polynomial is zero.
"""
iszero(p::PolynomialSparse)::Bool = p.terms == [Term(0,0)]

#isone(p::Polynomial)::Bool = p.terms == [Term()]
#################################################################
# Transformation of the polynomial to create another polynomial #
#################################################################

"""
The negative of a polynomial.
"""
-(p::PolynomialSparse) = PolynomialSparse(map((pt)->-pt, p.terms))

"""
Create a new polynomial which is the derivative of the polynomial.
"""
function derivative(p::PolynomialSparse)::PolynomialSparse 
    der_p = PolynomialSparse()
    for term in p
        der_term = derivative(term)
        !iszero(der_term) && push!(der_p,der_term)
    end
    return trim!(der_p)
end

"""
The prim part (multiply a polynomial by the inverse of its content). 
"""
############################################################################################################
prim_part(p::PolynomialSparse) = p ÷ content(p)

"""
A square free polynomial.
"""
############################################################################################################
square_free(p::Polynomial, prime::Int)::Polynomial = (p ÷ gcd(p,derivative(p),prime))(prime)

#################################
# Queries about two polynomials #
#################################

"""
Check if two polynomials are the same
"""
==(p1::PolynomialSparse, p2::PolynomialSparse)::Bool = p1.terms == p2.terms


"""
Check if a polynomial is equal to 0.
"""
#Note that in principle there is a problem here. E.g The polynomial 3 will return true to equalling the integer 2.
==(p::PolynomialSparse, n::T) where T <: Real = iszero(p) == iszero(n)

##################################################################
# Operations with two objects where at least one is a polynomial #
##################################################################

"""
Subtraction of two polynomials.
"""
#-(p1::PolynomialSparse, p2::PolynomialSparse)::PolynomialSparse = p1 + (-p2)
function -(p::PolynomialSparse, t::Term)::PolynomialSparse
    pp = deepcopy(p)
    cnt = 0
    for terms in pp
        cnt += 1
        if terms.degree == t.degree
            pp.terms[cnt] += (-t)
        else
            pp.terms[cnt] = terms
        end
    end

    if t.degree ∉ degrees(p)
        push!(pp, (-t))
    end
    return pp
end

function -(p1::PolynomialSparse, p2::PolynomialSparse)::PolynomialSparse
    p = deepcopy(p1)
    for t in p2
        p = p - t
    end
    return p
end

"""
Multiplication of polynomial and term.
"""
*(t::Term, p1::PolynomialSparse)::PolynomialSparse = iszero(t) ? PolynomialSparse() : PolynomialSparse(map((pt)->t*pt, p1.terms))
*(p1::PolynomialSparse, t::Term)::PolynomialSparse = t*p1

"""
Multiplication of polynomial and an integer.
"""
*(n::Int, p::PolynomialSparse)::PolynomialSparse = p*Term(n,0)
*(p::PolynomialSparse, n::Int)::PolynomialSparse = n*p

"""
Integer division of a polynomial by an integer.

Warning this may not make sense if n does not divide all the coefficients of p.
"""
÷(p::PolynomialSparse, n::Int) = (prime)->PolynomialSparse(map((pt)->((pt ÷ n)(prime)), p.terms))

"""
Take the mod of a polynomial with an integer.
"""
function mod(f::PolynomialSparse, p::Int)::PolynomialSparse
    f_out = deepcopy(f)
    for i in 1:length(f_out.terms)
        f_out.terms[i] = mod(f_out.terms[i], p)
    end
    return trim!(f_out)
        
    # p_out = Polynomial()
    # for t in f
    #     new_term = mod(t, p)
    #     @show new_term
    #     push!(p_out, new_term)
    # end
    # return p_out
end

"""
Power of a polynomial mod prime.
"""
function pow_mod(p::PolynomialSparse, n::Int, prime::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
        out = mod(out, prime)
    end
    return out
end

