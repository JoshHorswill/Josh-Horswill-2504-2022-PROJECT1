struct PolynomialModP
    
    p::PolynomialSparse
    prime::Int
    
    #terms::Vector{Term}   
    
    #Inner constructor of 0 polynomial
    #PolynomialSparse() = new([Term(0,0)])
    PolynomialModP() = new(PolynomialSparse(zero(Term)), 101)
    
    #Inner constructor of polynomial based on arbitrary list of terms
    function PolynomialModP(vt::PolynomialSparse, prime::Int)

        #Filter the vector so that there is not more than a single zero term
        #vt = filter((t)->!iszero(t), vt)
        if isempty(vt)
            vt = [zero(Term)]
        end

        # Take count of elements in vector of terms
        cnt = 0
        for t in vt
            cnt += 1
        end

        terms = [zero(Term) for i in 0:cnt-1] #Array of zeroes to the size of the count

        cnter = 0
        inner_ct = 0
        for t in vt
            inner_ct += 1
            if t.coeff !=0 && t.degree != 0
                cnter += 1
                #terms[cnter] = mod(t, prime) #+1 accounts for 1-indexing
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

        modp_sorted = mod(PolynomialSparse(terms), prime)
        return new(modp_sorted, prime)
    end
end


"""
Construct a polynomial with a single term.
"""
PolynomialModP(t::Term, prime::Int) = PolynomialModP(PolynomialSparse(t), prime)

"""
Construct a polynomial of the form x.
"""
x_poly_mod() = PolynomialModP(Term(1,1), 3)

"""
Creates the zero polynomial.
"""
zero(::Type{PolynomialModP})::PolynomialModP = PolynomialModP(PolynomialSparse(Term(0,0)), 3)

"""
Creates the unit polynomial.
"""
one(::Type{PolynomialModP})::PolynomialModP = PolynomialModP(PolynomialSparse(Term(1,1)), 3)
one(p::PolynomialModP)::PolynomialModP = PolynomialModP(PolynomialSparse(Term(1,1)), 3)


"""
Generates a random PolynomialModP
"""
function rand(::Type{PolynomialModP})
    primes = [1, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101]
    n = length(primes)
    x = rand(PolynomialSparse)
    idx = rand(1:n)
    pri = primes[idx]
    return PolynomialModP(x, 101)
end


"""
Show a polynomial.
"""
function show(io::IO, f::PolynomialModP) 
    print("( ")
    if iszero(f.p)
        print(io,"0")
    else
    
        n = length(f.p.terms)
        
        for (i,t) in enumerate(f.p.terms)
            if !iszero(t)
                if t.degree == 0 && t.coeff < 0
                    print(io, i != 1 ? " " : "", t.coeff, " ")
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
        println(" ) mod ", f.prime)
    end
end


function trim!(f::PolynomialModP)::PolynomialModP
    
    i = length(f.p.terms)
    while i > 1
        if iszero(f.p.terms[i])
            pop!(f.p.terms)
        else
            break
        end
        i -= 1
    end
    return f
end
##############################################
# Iteration over the terms of the polynomial mod p #
##############################################

"""
Allows to do iteration over the non-zero terms of the polynomial. This implements the iteration interface.
"""
iterate(f::PolynomialModP, state=1) = iterate(f.p.terms, state)

##############################
# Queries about a polynomial mod p#
##############################

"""
The number of terms of the polynomial.
"""
length(f::PolynomialModP) = length(f.p.terms) 

"""
The leading term of the polynomial.
"""
leading(f::PolynomialModP)::Term = isempty(f.p.terms) ? zero(Term) : first(f.p.terms) 

"""
Returns the coefficients of the polynomial.
"""
coeffs(f::PolynomialModP)::Vector{Int} = [t.coeff for t in f.p]

"""
The degree of the polynomial.
"""
degree(f::PolynomialModP)::Int = leading(f.p).degree 
degrees(f::PolynomialModP)::Vector{Int} = [t.degree for t in f.p]

"""
Push a new term into the polynomial.
"""
#Note that ideally this would throw and error if pushing another term of degree that is already in the polynomial
function push!(f::PolynomialModP, t::Term) 
    #if t.degree <= degree(p)
    #    p.terms[1] = t
    #else
    #    append!(p.terms, zeros(Term, t.degree - degree(p)-1))
    push!(f.p.terms, t)
    #end
    #return p        
end

"""
Pop the leading term out of the polynomial. When polynomial is 0, keep popping out 0.
"""
function pop!(f::PolynomialModP)::Term 
    popped_term = pop!(f.p.terms) #last element popped is leading coefficient

    while !isempty(f.p.terms) && iszero(last(f.p.terms))
        pop!(f.p.terms)
    end

    if isempty(f.p.terms)
        push!(f.p.terms, zero(Term))
    end

    return popped_term
end

"""
Check if the polynomial is zero.
"""
iszero(f::PolynomialModP)::Bool = f.p.terms == [Term(0,0)]

"""
The negative of a polynomial. #############################################
"""
-(f::PolynomialModP) = PolynomialModP(map((pt)->-pt, f.p.terms), f.prime)

"""
Create a new polynomial which is the derivative of the polynomial.
"""
==(f1::PolynomialModP, f2::PolynomialModP)::Bool = f1.p.terms == f2.p.terms
==(f::PolynomialModP, n::T) where T <: Real = iszero(f.p) == iszero(n)


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
    p = deepcopy(f1)
    for t in f2.p
        p += t
    end
    return p
end

+(p::PolynomialModP, n::Int) = p + Term(n,0)
+(n::Int, p::PolynomialModP) = p + Term(n,0)
-(p::PolynomialModP, t::Term) = p + (-t)

x = rand(PolynomialModP)
x.prime
t = Term(27,1)
-(x,t)

*(t::Term, p1::PolynomialModP)::PolynomialModP = iszero(t) ? PolynomialModP() : PolynomialModP(PolynomialSparse(map((pt)->t*pt, p1.p.terms)), p1.prime)
*(p1::PolynomialModP, t::Term)::PolynomialModP = t*p1
*(n::Int, p::PolynomialModP)::PolynomialModP = p*Term(n,0)
*(p::PolynomialModP, n::Int)::PolynomialModP = n*p

function *(p1::PolynomialModP, p2::PolynomialModP)::PolynomialModP
    p_out = PolynomialSparse()
    for t in p1.p
        new_summand = (t * p2.p)
        p_out += new_summand
    end
    return PolynomialModP(p_out, p1.prime)
end


function -(p::PolynomialModP, t::Term)::PolynomialModP
    pp = deepcopy(p.p)
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

    result = PolynomialModP(pp, p.prime)
    return result
end



function -(p1::PolynomialModP, p2::PolynomialModP)::PolynomialModP
    @assert p1.prime == p2.prime
    p = deepcopy(p1.p)
    for t in p2.p
        p = p - t
    end
    result = p
    return PolynomialModP(result, p1.prime)
end

PolynomialModP(PolynomialSparse([Term(1,1), Term(3,3)]), 101) - Term(2,1)

function mod(f::PolynomialModP, p::Int)::PolynomialModP
    f_out = deepcopy(f.p)
    for i in 1:length(f_out.terms)
        f_out.terms[i] = mod(f_out.terms[i], p)
    end
    return PolynomialModP(f_out, p)
        
    # p_out = Polynomial()
    # for t in f
    #     new_term = mod(t, p)
    #     @show new_term
    #     push!(p_out, new_term)
    # end
    # return p_out
end



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
Power of a polynomial mod prime.
"""
function pow_mod(p::PolynomialModP, n::Int)
    n < 0 && error("No negative power")
    out = one(p.p)
    for _ in 1:n
        out *= p.p
        out = mod(out, p.prime)
    end
    return out
end
