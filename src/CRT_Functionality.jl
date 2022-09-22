##### CHINESE REMAINDER THEOREM #####

### GET INDEX FUNCTION
function my_getindex(idx::Int64, v1::PolynomialSparse)::Int64
    v11 = degrees(v1)
    xarr = findall(x->x==idx, v11)
    return xarr[1]
end


### CRT OF TWO ARRAYS OF NUMBERS (SIZE N = 2)
function iCRT(um::Vector{Int64}, p::Vector{Int64})
    n = length(um)
    v = Vector{Int}(undef,n)

    v[1] = um[1]
    v[2] = (um[2] - v[1]) * Base.invmod(p[1], p[2]) % p[2]
    @show v #intermediate 
    u = v[1] + v[2]*p[1]

    return u
end

### CRT OF TWO GIVEN POLYNOMIALS
function CRT_poly(p1::PolynomialSparse, p2::PolynomialSparse, n::Int64, m::Int64)
    ak, bk = 0, 0
    c = PolynomialSparse(Term(0,0))
    lenmax = max(degree(p1), degree(p2))
    for k in 1:lenmax
        if k ∉ degrees(p1)
            ak = 0
        else
            spot = my_getindex(k, p1)
            ak = coeffs(p1)[spot]
        end

        if k ∉ degrees(p2)
            bk = 0
        else
            spot = my_getindex(k, p2)
            bk = coeffs(p2)[spot]
        end

        ck = iCRT([ak, bk], [n,m])
        c = c + Term(ck, k)
    end

    return c

end

### PRODUCT OF TWO GIVEN POLYNOMIALS USING CRT
using Primes
function CRT_prod(p1::PolynomialSparse, p2::PolynomialSparse)
    height_a = maximum(coeffs(p1))
    height_b = maximum(coeffs(p2))
    
    bound = 2 * height_a * height_b * min(degree(p1)+1, degree(p2)+1)
    pri = 3
    M = pri

    fir_prod = p1*p2

    C = mod(fir_prod, M)


    while M < bound
        pri = nextprime(pri+1)
        Cdash = mod(p1*p2, pri)
        C = CRT_poly(C, Cdash, M, pri)
        M = M*pri
    end

    return C
end


### CRT FULL IMPLEMENTATION FOR BINT

function my_getindex(idx::Int64, v1::PolynomialSparseBInt)::Int64
    v11 = degrees(v1)
    xarr = findall(x->x==idx, v11)
    return xarr[1]
end

function iCRT_BI(um::Vector{BigInt}, p::Vector{BigInt})
    n = length(um)
    v = Vector{BigInt}(undef,n)

    v[1] = um[1]
    v[2] = (um[2] - v[1]) * Base.invmod(p[1], p[2]) % p[2]
    u = v[1] + v[2]*p[1]

    return u
end

function CRT_polyBI(p1::PolynomialSparseBInt, p2::PolynomialSparseBInt, n::BigInt, m::BigInt)
    ak, bk = 0, 0
    c = PolynomialSparseBInt(TermBI(big"0",0))
    lenmax = max(degree(p1), degree(p2))
    for k in 1:lenmax
        if k ∉ degrees(p1)
            ak = 0
        else
            spot = my_getindex(k, p1)
            ak = coeffs(p1)[spot]
        end

        if k ∉ degrees(p2)
            bk = 0
        else
            spot = my_getindex(k, p2)
            bk = coeffs(p2)[spot]
        end

        ck = iCRT_BI([BigInt(ak), BigInt(bk)], [n,m])
        c = c + TermBI(BigInt(ck), k)
    end
    return c
end


function CRT_prodBI(p1::PolynomialSparseBInt, p2::PolynomialSparseBInt)
    height_a = maximum(coeffs(p1))
    height_b = maximum(coeffs(p2))
    
    bound = 2 * height_a * height_b * min(degree(p1)+1, degree(p2)+1)
    pri = big"3"
    M = pri

    C = mod(p1, BigInt(M))*mod(p2, BigInt(M))

    while M < bound
        pri = nextprime(pri+1)
        Cdash = mod(p1, BigInt(pri))*mod(p2, BigInt(pri))
        C = CRT_polyBI(C, Cdash, M, pri)
        M = M*pri
    end

    return C
end


