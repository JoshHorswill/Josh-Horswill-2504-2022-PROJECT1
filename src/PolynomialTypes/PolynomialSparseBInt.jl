### PolynomialSparseBInt: development of type and functionality

struct PolynomialSparseBInt
    
    terms::Vector{TermBI}   
    
    #Inner constructor of 0 polynomial
    PolynomialSparseBInt() = new([zero(TermBI)])

    #Inner constructor of polynomial based on arbitrary list of terms
    function PolynomialSparseBInt(vt::Vector{TermBI})

        #Filter and sort zero TermBI
        vt = filter((t)->!iszero(t), vt)
        if isempty(vt)
            vt = [zero(TermBI)]
        end

        # Take count of elements in vector of terms
        cnt = 0
        for t in vt
            cnt += 1
        end

        terms = [zero(TermBI) for i in 0:cnt-1] #Array of zeroes to the size of the count

        # For each term, if not a null item, commit to terms array and increment index
        cnter = 0
        for t in vt
            if t.coeff !=0 && t.degree != 0 
                cnter += 1
                terms[cnter] = t
            elseif t.coeff != 0 && t.degree == 0
                cnter += 1
                terms[cnter] = t
            end
        end

        # Sort by degree and return
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

function trim!(p::PolynomialSparseBInt)::PolynomialSparseBInt
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

### SINGLE TERM CONSTRUCTOR
PolynomialSparseBInt(t::TermBI) = PolynomialSparseBInt([t])

### CYCLOTONIC CONSTRUCTOR
cyclotonic_polynomial(p::Int) = PolynomialSparseBInt([TermBI(big"1",p), TermBI(big"-1",0)])

### LINEAR MONIC CONSTRUCTOR
linear_monic_polynomial(n::Int) = PolynomialSparseBInt([TermBI(big"1",1), TermBI(big"-n",0)])

### SINGLE VARIABLE DEFINITION
x_poly_BI() = PolynomialSparseBInt(TermBI(big"1",1))

### ZERO TERMBI
zero(::Type{PolynomialSparseBInt})::PolynomialSparseBInt = PolynomialSparseBInt(TermBI(big"0",0))

### UNIT POLYNOMIAL OF TERMBI
one(::Type{PolynomialSparseBInt})::PolynomialSparseBInt = PolynomialSparseBInt(TermBI(big"1",0))
one(p::PolynomialSparseBInt)::PolynomialSparseBInt = PolynomialSparseBInt(TermBI(big"1",0))

### RANDOM BIGINT POLYNOMIAL GENERATOR
function rand(::Type{PolynomialSparseBInt} ; 
    degree::Int = -1, 
    terms::Int = -1, 
    max_coeff::BigInt = big"10", 
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
        p = PolynomialSparseBInt( [TermBI(coeffs[i],degrees[i]) for i in 1:length(degrees)] )
        condition(p) && return p
    end
end


### SHOW POLYNOMIALSPARSEBINT W/ PRETTY PRINTING
function show(io::IO, p::PolynomialSparseBInt) 
    if iszero(p)
        print(io,"0")
    else
        n = length(p.terms)
        for (i,t) in enumerate(p.terms)
            if !iszero(t)
                if t.degree == 0 && t.coeff < 0
                    print(io, i != 1 ? " " : "", t.coeff, " ")
                elseif t.degree == 1 && t.coeff > 0
                    print(io, i != 1 ? " + " : "", t.coeff, "x")
                elseif t.degree == 1 && t.coeff < 0
                    print(io, i != 1 ? " " : "", t.coeff, "x")
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


### IS POLYNOMIAL 0: Used a different function for this type due to BigInt usage
function iszero2(p::PolynomialSparseBInt)::Bool
    ret = false
    for term in p
        if term.coeff == big"0" & term.degree == 0
            ret = true
        end
    end
    return ret
end

# Original: kept for comparison
iszero(p::PolynomialSparseBInt)::Bool = p.terms == [TermBI(big"0",0)]

### ITERATE POLYNOMIAL TERMS
iterate(p::PolynomialSparseBInt, state=1) = iterate(p.terms, state)

### LENGTH OF POLYNOMIAL
length(p::PolynomialSparseBInt) = length(p.terms) 

### LEADING TERM OF POLYNOMIAL
function leading(p::PolynomialSparseBInt)::TermBI
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
    return TermBI(co, de)
end

### COEFFICIENTS AND DEGREES OF POLYNOMAIL
coeffs(p::PolynomialSparseBInt)::Vector{Int} = [t.coeff for t in p]
degrees(p::PolynomialSparseBInt)::Vector{Int} = [t.degree for t in p]

function degree(p::PolynomialSparseBInt)::Int
    max_deg = -Inf
    for terms in p
        if terms.degree > max_deg
            max_deg = terms.degree
        end
    end
    return max_deg
end

### EUCLIDEAN CONTENT OF POYLNOMIAL
content(p::PolynomialSparseBInt)::Int = euclid_alg(coeffs(p))

### EVALUATOR POLYNOMIAL AT GIVEN POINT
evaluate(f::PolynomialSparseBInt, x::T) where T <: Number = sum(evaluate(t,x) for t in f)

### PUSH TERM INTO POLYNOMIALSPARSEBINT
function push!(p::PolynomialSparseBInt, t::TermBI) 
    push!(p.terms, t)
end

### POP TERM FROM POLYNOMIALSPARSEBINT
# Not used
function pop!(p::PolynomialSparseBInt)::TermBI
    popped_term = pop!(p.terms) #last element popped is leading coefficient

    while !isempty(p.terms) && iszero(last(p.terms))
        pop!(p.terms)
    end

    if isempty(p.terms)
        push!(p.terms, zero(TermBI))
    end

    return popped_term
end

### NEGATION OF POLYNOMIALSPARSEBINT
-(p::PolynomialSparseBInt) = PolynomialSparseBInt(map((pt)->-pt, p.terms))

### DERIVATIVE OF POLYNOMIALSPARSEBINT
function derivative(p::PolynomialSparseBInt)::PolynomialSparseBInt 
    der_p = PolynomialSparseBInt()
    for term in p
        der_term = derivative(term)
        !iszero(der_term) && push!(der_p,der_term)
    end
    return trim!(der_p)
end

### MODULO OF POLYNOMIALSPARSEBINT
function mod(f::PolynomialSparseBInt, p::BigInt)::PolynomialSparseBInt
    f_out = deepcopy(f)
    for i in 1:length(f_out.terms)
        f_out.terms[i] = mod(f_out.terms[i], p)
    end
    return trim!(f_out)
end

### SYMMETRIC MODULO OF POLYNOMIALSPARSEBINT
function smod(f::PolynomialSparseBInt, p::BigInt)::PolynomialSparseBInt
    f_out = deepcopy(f)
    for i in 1:length(f_out.terms)
        if mod(i, p) <= p//2
            f_out.terms[i] = mod(f_out.terms[i], p)
        else
            f_out.terms[i] = mod(f_out.terms[i], p) - p
        end
    end

    return trim!(f_out)

end


### PRIME PART OF POLYNOMIAL
prim_part(p::PolynomialSparseBInt) = p ÷ content(p)


### EQUATE OR CHECK 0 FOR POLYNOMIALS
==(p1::PolynomialSparseBInt, p2::PolynomialSparseBInt)::Bool = p1.terms == p2.terms
==(p::PolynomialSparseBInt, n::T) where T <: Real = iszero(p) == iszero(n)


### POLYNOMIALSPARSEBINT SUBTRACTION
# Again, much rather its own function than using addition
function -(p::PolynomialSparseBInt, t::TermBI)::PolynomialSparseBInt
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

function -(p1::PolynomialSparseBInt, p2::PolynomialSparseBInt)::PolynomialSparseBInt
    p = deepcopy(p1)
    for t in p2
        p = p - t
    end
    return p
end


### TERMBI AND BIGINT MULTIPLICATION FOR POLYNOMIALSPARSEBINT
*(t::TermBI, p1::PolynomialSparseBInt)::PolynomialSparseBInt = iszero(t) ? PolynomialSparseBInt() : PolynomialSparseBInt(map((pt)->t*pt, p1.terms))
*(p1::PolynomialSparseBInt, t::TermBI)::PolynomialSparseBInt = t*p1
*(n::BigInt, p::PolynomialSparseBInt)::PolynomialSparseBInt = p*TermBI(n,0)
*(p::PolynomialSparseBInt, n::BigInt)::PolynomialSparseBInt = n*p


### DIVISION INDICATOR
÷(p::PolynomialSparseBInt, n::BigInt) = (prime)->PolynomialSparse(map((pt)->((pt ÷ n)(prime)), p.terms))


### POWER MOD OF POLYNOMIALSPARSEBINT
function pow_mod(p::PolynomialSparseBInt, n::BigInt, prime::BigInt)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
        out = mod(out, prime)
    end
    return out
end

### REDEFINITION W/ REPEATED SQUARING
function pow_mod2(p::PolynomialSparseBInt, n::Int, prime::BigInt)
    n < 0 && error("No negative power")
    out = one(p)
    retur = repsqr_exp(p, n)
    #modret = mod(retur, prime)
    modret = retur

    return modret
end

