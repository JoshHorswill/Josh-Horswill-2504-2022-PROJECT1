"""
A Polynomial Sparse Type: Conservative Dense Approach
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

        bubble_poly_sort!(terms)

        return new(terms)
    end
end

### DEGREE SORT FUNCTION
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


### TRIM IN PLACE FUNCTION
#This function maintains the invariant of the Polynomial type so that there are no zero terms beyond the highest
#non-zero term.

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



### POLYNOMIALSPARSE TERM CONSTRUCTOR
PolynomialSparse(t::Term) = PolynomialSparse([t])

### CYCLOTONIC SPARSE 
cyclotonic_polynomial(p::Int) = PolynomialSparse([Term(1,p), Term(-1,0)])

### LINEAR MONIC SPARSE
linear_monic_polynomial(n::Int) = PolynomialSparse([Term(1,1), Term(-n,0)])

### POLYNOMIAL SPARSE X DEFINITION
x_poly_sp() = PolynomialSparse(Term(1,1))

### ZERO POLYNOMIALSPARSE DEFINITION
zero(::Type{PolynomialSparse})::PolynomialSparse = PolynomialSparse(Term(0,0))

### UNIT POLYNOMIALSPARSE DEFINITION: CREATES A 1
one(::Type{PolynomialSparse})::PolynomialSparse = PolynomialSparse(Term(1,0))
one(p::PolynomialSparse)::PolynomialSparse = PolynomialSparse(Term(1,0))


### RANDOM POLYNOMIALSPARSE GENERATOR
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

### SHOW A POLYNOMIAL ("PRETTY PRINTING")
function show(io::IO, p::PolynomialSparse) 
    if iszero(p)
        print(io,"0")
    else
        n = length(p.terms)
        for (i,t) in enumerate(p.terms)
            if !iszero(t)
                if t.degree == 0 && t.coeff < 0
                    print(io, i != 1 ? " " : "", t.coeff, "x")
                elseif t.degree == 1 && t.coeff > 0
                    print(io, i != 1 ? " + " : "", t.coeff, "x")
                elseif t.degree == 1 && t.coeff < 0
                    print(io, i != 1 ? " " : "", t.coeff, "x")
                elseif t.degree == 0 && t.coeff > 0
                    print(io, i != 1 ? " + " : "", t.coeff, "", "")
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

### ITERATE OVER SPARSE TERMS IN POLYNOMIAL
iterate(p::PolynomialSparse, state=1) = iterate(p.terms, state)




##############################
# Queries about a polynomial #
##############################

### LENGTH OF SPARSE TERMS IN POLYNOMIAL
length(p::PolynomialSparse) = length(p.terms) 

### LEADING TERM OF POLYNOMIAL (RE-WORKED)
function leading(p::PolynomialSparse)::Term
    lead = -Inf
    co, de = 0, 0
    for terms in p
        if terms.degree > lead
            lead = terms.degree
            co = terms.coeff
            de = terms.degree
        end
    end
    return Term(co, de)
end

### COEFFICIENTS & DEGREES IN POLYNOMIALSPARSE COLLECTION
coeffs(p::PolynomialSparse)::Vector{Int} = [t.coeff for t in p]
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


### EUCLIDEAN CONTENT OF POLYNOMIALSPARSE
content(p::PolynomialSparse)::Int = euclid_alg(coeffs(p))

### EVALUATOR AT A GIVEN POINT
evaluate(f::PolynomialSparse, x::T) where T <: Number = sum(evaluate(t,x) for t in f)

################################
# Pushing and popping of terms #
################################

### PUSH FUNCTION
# Using the base function and manipulating when needed worked better
function push!(p::PolynomialSparse, t::Term) 
    push!(p.terms, t)
end

### POP FUNCTION
# Not used
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

### ISZERO CHECK
iszero(p::PolynomialSparse)::Bool = p.terms == [Term(0,0)]


#################################################################
# Transformation of the polynomial to create another polynomial #
#################################################################

### NEGATING POLYNOMIAL
-(p::PolynomialSparse) = PolynomialSparse(map((pt)->-pt, p.terms))

### POLYNOMIALSPARSE DERIVATIVE
function derivative(p::PolynomialSparse)::PolynomialSparse 
    der_p = PolynomialSparse()
    for term in p
        der_term = derivative(term)
        !iszero(der_term) && push!(der_p,der_term)
    end
    return trim!(der_p)
end

### PRIME PART OF POLYNOMIALSPARSE
prim_part(p::PolynomialSparse) = p ÷ content(p)

### SQUARE FREE POLYNOMIALSPARSE
square_free(p::PolynomialSparse, prime::Int)::PolynomialSparse = (p ÷ gcd(p,derivative(p),prime))(prime)



#################################
# Queries about two polynomials #
#################################

### EQUALITY CHECK OF TWO POLYNOMIALSPARSES
==(p1::PolynomialSparse, p2::PolynomialSparse)::Bool = p1.terms == p2.terms

### IS POLYNOMIALSPARSE 0 FUNCTION
# not used
==(p::PolynomialSparse, n::T) where T <: Real = iszero(p) == iszero(n)

##################################################################
# Operations with two objects where at least one is a polynomial #
##################################################################

### SUBTRACTION OF TWO POLYNOMIAL SPARSES
# Preferred using its own function 
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
    kar = p2.terms
    for t in p2
        p = -(p, t)
    end
    return p
end

function -(p1::PolynomialSparse, p2::PolynomialSparse)::PolynomialSparse
    p = deepcopy(p1)
    if p1 == p2
        return PolynomialSparse(Term(0,0))
    else
        for t in p2
            p = PolynomialSparse((p - t).terms)
        end
        return p
    end
end


### TERM MULTIPLICATION
*(t::Term, p1::PolynomialSparse)::PolynomialSparse = iszero(t) ? PolynomialSparse() : PolynomialSparse(map((pt)->t*pt, p1.terms))
*(p1::PolynomialSparse, t::Term)::PolynomialSparse = t*p1

### INTEGER MULTIPLICATION
*(n::Int, p::PolynomialSparse)::PolynomialSparse = p*Term(n,0)
*(p::PolynomialSparse, n::Int)::PolynomialSparse = n*p

### DIVISION OPERATOR
÷(p::PolynomialSparse, n::Int) = (prime)->PolynomialSparse(map((pt)->((pt ÷ n)(prime)), p.terms))

### MOD OPERATOR
function mod(f::PolynomialSparse, p::Int)::PolynomialSparse
    f_out = deepcopy(f)
    for i in 1:length(f_out.terms)
        f_out.terms[i] = mod(f_out.terms[i], p)
    end
    return trim!(f_out)
end

### POWER MOD OPERATOR
function pow_mod(p::PolynomialSparse, n::Int, prime::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
        out = mod(out, prime)
    end
    return out
end

