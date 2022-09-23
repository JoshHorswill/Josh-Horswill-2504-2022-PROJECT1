struct PolynomialModP
    # POLYNOMIALMODP: Takes in PolynomialSparse and Prime, performs sparse functionality with respect to modulo
    
    # Intake
    p::PolynomialSparse
    prime::Int
    
    # Zero case, prime is given to be 101
    PolynomialModP() = new(PolynomialSparse(zero(Term)), 101)
    
    #Inner constructor of polynomial based on arbitrary list of terms
    function PolynomialModP(vt::PolynomialSparse, prime::Int)

        # Set vector of terms if empty
        if isempty(vt)
            vt = [zero(Term)]
        end

        # Take count of elements in vector of terms
        cnt = 0
        for t in vt
            cnt += 1
        end

        terms = [zero(Term) for i in 0:cnt-1] #Array of zeroes to the size of the count

        # Increment and set polynomial terms
        cnter = 0
        inner_ct = 0
        for t in vt
            inner_ct += 1
            if t.coeff !=0 && t.degree != 0
                cnter += 1
                terms[cnter] = t
            elseif t.coeff != 0 && t.degree == 0
                cnter += 1
                terms[cnter] = t
            end
        end

        # Sort w.r.t descending degree
        bubble_poly_sort!(terms)

        # Set and return polynomial, prime
        modp_sorted = PolynomialSparse(terms)
        return new(modp_sorted, prime)
    end
end

# Sort for degrees function 
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


### SINGLE TERM MOD P
PolynomialModP(t::Term, prime::Int) = PolynomialModP(PolynomialSparse(t), prime)

### X FORM MOD P
x_poly_mod() = PolynomialModP(Term(1,1), 3)

### ZERO TERM MOD P
zero(::Type{PolynomialModP})::PolynomialModP = PolynomialModP(PolynomialSparse(Term(0,0)), 3)

### UNIT TERM MOD P
one(::Type{PolynomialModP})::PolynomialModP = PolynomialModP(PolynomialSparse(Term(1,1)), 3)
one(p::PolynomialModP)::PolynomialModP = PolynomialModP(PolynomialSparse(Term(1,1)), 3)


### RANDOM FUNCTION MOD P
# Uses PolynomialSparse random, resemove pri variable and insert idx to output to make "random",
# the random prime is set to 5 currently to allow consistent operation.
function rand(::Type{PolynomialModP})
    primes = [1, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101]
    n = length(primes)
    x = rand(PolynomialSparse)
    idx = rand(1:n)
    pri = primes[3]
    return PolynomialModP(x, pri)
end


### SHOW MOD P 
function show(io::IO, f::PolynomialModP) 
    print("( ")
    if iszero(f.p)
        print(io,"0")
    else
    
        n = length(f.p.terms)
        pmod = mod(f.p, f.prime)

        for (i,t) in enumerate(pmod)
            if !iszero(t)
                if t.degree == 0 && t.coeff < 0
                    print(io, i != 1 ? " " : "", t.coeff, " ")
                elseif t.degree == 1 && t.coeff > 0
                    print(io, i != 1 ? " + " : "", t.coeff, "x")
                elseif t.degree == 1 && t.coeff < 0
                    print(io, i != 1 ? " " : "", t.coeff, "x")
                elseif t.degree == 1 && t.coeff > 0
                    print(io, i != 1 ? " + " : "", t.coeff, "x")
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




### TRIM AND SET MODP
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

### ITERATE THROUGH TERMS 
iterate(f::PolynomialModP, state=1) = iterate(f.p.terms, state)

##############################
# Queries about a polynomial mod p#
##############################

### LENGTH OF MODP POLYNOMIAL
length(f::PolynomialModP) = length(f.p.terms) 

### LEADING TERM OF MODP POLYNOMIAL
leading(f::PolynomialModP)::Term = isempty(f.p.terms) ? zero(Term) : first(f.p.terms) 

### COEFFICIENTS OF MODP POLYNOMIAL 
coeffs(f::PolynomialModP)::Vector{Int} = [t.coeff for t in f.p]

### DEGREE & DEGREES OF MODP POLYNOMIAL
degree(f::PolynomialModP)::Int = leading(f.p).degree 
degrees(f::PolynomialModP)::Vector{Int} = [t.degree for t in f.p]

### PUSH TERM INTO MODP POLYNOMIAL
# Use sparringly
function push!(f::PolynomialModP, t::Term) 
    push!(f.p.terms, t)   
end


### POP TERM FROM POLYNOMIAL MODP
# Not used
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

### IS ZERO CHECK POLYNOMIAL MODP
iszero(f::PolynomialModP)::Bool = f.p.terms == [Term(0,0)]

### NEGATIVE POLYNOMIAL MODP
-(f::PolynomialModP) = PolynomialModP(map((pt)->-pt, f.p.terms), f.prime)

### EQUATE POLYNOMIAL MODP
==(f1::PolynomialModP, f2::PolynomialModP)::Bool = f1.p.terms == f2.p.terms
==(f::PolynomialModP, n::T) where T <: Real = iszero(f.p) == iszero(n)

### POLYNOMIAL MODP SUBTRACTION:
# Prefered to use a function of its own instead of compliment addition function
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

    result = PolynomialModP(mod(pp, p.prime), p.prime)
    return result
end


### ( POLYNOMIAL ON POLYNOMIAL ) MOD P FUNCTION
function -(p1::PolynomialModP, p2::PolynomialModP)::PolynomialModP
    @assert p1.prime == p2.prime
    p = deepcopy(p1.p)
    for t in p2.p
        p = p - t
    end
    result = p
    return PolynomialModP(mod(result, p1.prime), p1.prime)
end

### POLYNOMIAL MOD P MOD FUNCTION
# Not used, PolynomialSparse Mod was all taht was needed
function mod(f::PolynomialModP, p::Int)::PolynomialModP
    f_out = deepcopy(f.p)
    for i in 1:length(f_out.terms)
        f_out.terms[i] = mod(f_out.terms[i], p)
    end
    return PolynomialModP(f_out, p)
end




"""
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
"""