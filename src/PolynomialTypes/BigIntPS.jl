struct PolynomialSparseBInt
    
    terms::Vector{TermBI}   
    
    #Inner constructor of 0 polynomial
    PolynomialSparseBInt() = new([zero(TermBI)])

    #Inner constructor of polynomial based on arbitrary list of terms
    function PolynomialSparseBInt(vt::Vector{TermBI})

        #Filter the vector so that there is not more than a single zero term
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

"""
Construct a polynomial with a single term.
"""
PolynomialSparseBInt(t::TermBI) = PolynomialSparseBInt([t])

"""
Construct a polynomial of the form x^p-x.
"""
cyclotonic_polynomial(p::Int) = PolynomialSparseBInt([TermBI(big"1",p), TermBI(big"-1",0)])


"""
Construct a polynomial of the form x-n.
"""
linear_monic_polynomial(n::Int) = PolynomialSparseBInt([TermBI(big"1",1), TermBI(big"-n",0)])

"""
Construct a polynomial of the form x.
"""
x_poly_BI() = PolynomialSparseBInt(TermBI(big"1",1))

"""
Creates the zero polynomial.
"""
#zero(::Type{PolynomialSparse})::PolynomialSparse = PolynomialSparse()
zero(::Type{PolynomialSparseBInt})::PolynomialSparseBInt = PolynomialSparseBInt(TermBI(big"0",0))

"""
Creates the unit polynomial.
"""
#one(::Type{PolynomialSparse})::PolynomialSparse = PolynomialSparse(one(Term))
#one(p::PolynomialSparse) = one(typeof(p))
one(::Type{PolynomialSparseBInt})::PolynomialSparseBInt = PolynomialSparseBInt(TermBI(big"1",0))
one(p::PolynomialSparseBInt)::PolynomialSparseBInt = PolynomialSparseBInt(TermBI(big"1",0))

"""
Generates a random polynomial.
"""

function rand(::Type{PolynomialSparseBInt} ; 
    degree::Int = -1, 
    terms::Int = -1, 
    max_coeff::BigInt = big"100", 
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

function show(io::IO, p::PolynomialSparseBInt) 
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


"""
Check if the polynomial is zero.
"""

function iszero2(p::PolynomialSparseBInt)::Bool
    ret = false
    for term in p
        if term.coeff == big"0" & term.degree == 0
            ret = true
        end
    end
    return ret
end
iszero(p::PolynomialSparseBInt)::Bool = p.terms == [TermBI(big"0",0)]
"""
Allows to do iteration over the non-zero terms of the polynomial. This implements the iteration interface.
"""
iterate(p::PolynomialSparseBInt, state=1) = iterate(p.terms, state)

"""
The number of terms of the polynomial.
"""
length(p::PolynomialSparseBInt) = length(p.terms) 

"""
The leading term of the polynomial.
"""
#leading(p::PolynomialSparse)::Term = isempty(p.terms) ? zero(Term) : first(p.terms) 
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

"""
Returns the coefficients of the polynomial.
"""
coeffs(p::PolynomialSparseBInt)::Vector{Int} = [t.coeff for t in p]


"""
The degree of the polynomial.
"""
#degree(p::PolynomialSparse)::Int = leading(p).degree 
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

"""
The content of the polynomial is the GCD of its coefficients.
"""
content(p::PolynomialSparseBInt)::Int = euclid_alg(coeffs(p))

"""
Evaluate the polynomial at a point `x`.
"""
evaluate(f::PolynomialSparseBInt, x::T) where T <: Number = sum(evaluate(t,x) for t in f)

function push!(p::PolynomialSparseBInt, t::TermBI) 
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

"""
The negative of a polynomial.
"""
-(p::PolynomialSparseBInt) = PolynomialSparseBInt(map((pt)->-pt, p.terms))

"""
Create a new polynomial which is the derivative of the polynomial.
"""
function derivative(p::PolynomialSparseBInt)::PolynomialSparseBInt 
    der_p = PolynomialSparseBInt()
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
prim_part(p::PolynomialSparseBInt) = p ÷ content(p)

"""
A square free polynomial.
"""
############################################################################################################
square_free(p::PolynomialSparseBInt, prime::Int)::PolynomialSparseBInt = (p ÷ gcd(p,derivative(p),prime))(prime)

"""
Check if two polynomials are the same
"""
==(p1::PolynomialSparseBInt, p2::PolynomialSparseBInt)::Bool = p1.terms == p2.terms


"""
Check if a polynomial is equal to 0.
"""
#Note that in principle there is a problem here. E.g The polynomial 3 will return true to equalling the integer 2.
==(p::PolynomialSparseBInt, n::T) where T <: Real = iszero(p) == iszero(n)


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


"""
Multiplication of polynomial and term.
"""
*(t::TermBI, p1::PolynomialSparseBInt)::PolynomialSparseBInt = iszero(t) ? PolynomialSparseBInt() : PolynomialSparseBInt(map((pt)->t*pt, p1.terms))
*(p1::PolynomialSparseBInt, t::TermBI)::PolynomialSparseBInt = t*p1


"""
Multiplication of polynomial and an integer.
"""
*(n::BigInt, p::PolynomialSparseBInt)::PolynomialSparseBInt = p*TermBI(n,0)
*(p::PolynomialSparseBInt, n::BigInt)::PolynomialSparseBInt = n*p


"""
Integer division of a polynomial by an integer.

Warning this may not make sense if n does not divide all the coefficients of p.
"""
÷(p::PolynomialSparseBInt, n::BigInt) = (prime)->PolynomialSparse(map((pt)->((pt ÷ n)(prime)), p.terms))

"""
Take the mod of a polynomial with an integer.
"""
function divide(num::PolynomialModP, den::PolynomialModP)
    @assert num.prime == den.prime
    if degree(den.p) > degree(num.p)
        return 0
    end

    num_val = num.p
    den_val = den.p
    println(typeof(num_val))
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



"""
Power of a polynomial mod prime.
"""
function pow_mod(p::PolynomialSparseBInt, n::BigInt, prime::BigInt)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
        out = mod(out, prime)
    end
    return out
end