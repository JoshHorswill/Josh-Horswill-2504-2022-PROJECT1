####################################
##### POLYNOMIALSPARSE TESTING #####
####################################
"""
Test product of polynomialsparse types
"""
function prod_test_poly_sparse(;N::Int = 10^3, N_prods::Int = 20, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        p2 = rand(PolynomialSparse)
        prod = p1*p2
        @assert p1.terms[1] * p2.terms[1] == prod.terms[2]
    end
    println("Polynomial Product - PASSED")
end


"""
Test derivative of polynomialssparse types
"""
function prod_derivative_test_poly_sparse(;N::Int = 10^2,  seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        p2 = rand(PolynomialSparse)
        p1d = derivative(p1)
        p2d = derivative(p2)
        @assert (p1d*p2 + p2d*p1).terms[2] == (derivative(p1*p2)).terms[2]
    end
    println("Polynomial Differentiation - PASSED")
end



"""
Test division of polynomialsparse modulo p.
"""
function division_test_poly(;prime::Int = 101, N::Int = 8, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        p2 = rand(PolynomialSparse)
        p_prod = p1*p2
        q, r = PolynomialSparse(), PolynomialSparse()
        try
            q, r = divide(p_prod, p2)(prime)
            if (q, r) == (nothing,nothing)
                println("Unlucky prime: $p1 is reduced to $(p1 % prime) modulo $prime")
                continue
            end
        catch e
            if typeof(e) == DivideError
                @assert mod(p2, prime) == 0
            else
                throw(e)
            end
        end
        @assert iszero( mod(q*p2+r - p_prod, prime) )
    end
    println("Polynomial Division - PASSED")
end
division_test_poly()


function ext_euclid_test_poly(;prime::Int=101, N::Int = 10^3, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        p2 = rand(PolynomialSparse)
        g, s, t = extended_euclid_alg(p1, p2, prime)
        @assert mod(s*p1 + t*p2 - g, prime) == 0
    end
    println("Polynomial Euclidean - PASSED")
end


#########################################
##### POLYNOMIALSPARSEBINT TESTING ######
#########################################

"""
Test product of polynomialsparse BigInt types
"""
function prod_test_poly_sparseBI(;N::Int = 10^3, N_prods::Int = 20, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparseBInt)
        p2 = rand(PolynomialSparseBInt)
        prod = p1*p2

        @assert coeffs(p1*p2)[end] == coeffs(prod)[end] && degrees(p1*p2)[end] == degrees(prod)[end]
    end
    println("Polynomial BI Product - PASSED")
end



"""
Test derivative of polynomials BigInt
"""
function prod_derivative_test_poly_sparseBI(;N::Int = 10^3,  seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparseBInt)
        p2 = rand(PolynomialSparseBInt)
        p1d = derivative(p1)
        p2d = derivative(p2)
        test11 = p1d*p2 + p2d*p1
        test12 = derivative(p1*p2)
        @assert coeffs(test11)[end] == coeffs(test12)[end] || coeffs(test11)[end] == coeffs(test12)[end-1]
    end
    println("Polynomial BI Differentiation - PASSED")
end



"""
Test division of polynomials BigInt modulo p.
"""
function division_test_polyBI(;prime::BigInt = big"101", N::Int = 10^3, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparseBInt)
        p2 = rand(PolynomialSparseBInt)
        p_prod = p1*p2
        q, r = PolynomialSparseBInt(), PolynomialSparseBInt()
        try
            q, r = divide(p_prod, p2)(prime)
            if (q, r) == (nothing,nothing)
                println("Unlucky prime: $p1 is reduced to $(p1 % prime) modulo $prime")
                continue
            end
        catch e
            if typeof(e) == DivideError
                @assert mod(p2, prime) == 0
            else
                throw(e)
            end
        end
        yy = mod(q*p2+r - p_prod, prime)
        @assert coeffs(yy)[1] == 0 
       
    end
    println("Polynomial BI Division - PASSED")
end


"""
Euclidean Algorithm testing for Polynomial BigInt
"""
function ext_euclid_test_polyBI(;prime::BigInt=big"101", N::Int = 5, seed::Int = 0)
    #Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparseBInt)
        p2 = rand(PolynomialSparseBInt)
        g, s, t = extended_euclid_alg(p1, p2, prime)
        xx = mod((s*p1 + t*p2) - g, prime)
        @assert degrees(xx)[1] == 0 && coeffs(xx)[1] == 0
    end
    println("Polynomial BI Euclidean - PASSED")
end


########################################
##### POLYNOMIALSPARSEMODP TESTING #####
########################################
"""
Test product of PolynomialModP types
"""
function prod_test_poly_modp(;N::Int = 10^3, N_prods::Int = 20, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialModP)
        p2 = rand(PolynomialModP)
        p1p = p1.p
        p2p = p2.p
        prodd = p1p*p2p
        @assert p1p.terms[end] * p2p.terms[end] == prodd.terms[end]
    end
    println("Polynomial ModP Product - PASSED")
end



"""
Test derivative of polynomialmodp
"""
function prod_derivative_test_poly_modp(;N::Int = 1,  seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialModP)
        p2 = rand(PolynomialModP)
        p1d = derivative(p1.p)
        p2d = derivative(p2.p)
        
        tes1 = p1d*p2.p + p2d*p1.p
        repack = PolynomialSparse([i for i in tes1])
        tes2 = derivative(p1.p*p2.p)
        @assert repack.terms[end] == tes2.terms[end]
    end
    println("Polynomial ModP Differentiation - PASSED")
end

### Note that primes are fixed for division adn euclidean testing seeing that ModP must be the same ###

"""
Test division of polynomialmodp
"""
function division_test_poly_modp(;prime::Int = 101, N::Int = 10, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialModP)
        p2 = rand(PolynomialModP)
        p1p = p1.p
        p2p = p2.p
        p_prod = p1.p*p2.p
        q, r = PolynomialModP(), PolynomialModP()
        try
            q, r = divide(p_prod, p2.p)(prime)
            if (q, r) == (nothing,nothing)
                println("Unlucky prime: $p1p is reduced to $(p1p % prime) modulo $prime")
                continue
            end
        catch e
            if typeof(e) == DivideError
                @assert mod(p2.p, prime) == 0
            else
                throw(e)
            end
        end
        @assert iszero( mod(q*p2.p+r - p_prod, prime) )
    end
    println("Polynomial ModP Division - PASSED")
end

"""
Euclidean Algorithm for PolynomialModP
"""
function ext_euclid_test_poly_modp(;prime::Int=101, N::Int = 10, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialModP)
        p2 = rand(PolynomialModP)
        p1p = p1.p
        p2p = p2.p
        g, s, t = extended_euclid_alg(p1p, p2p, prime)
        @assert mod(s*p1p + t*p2p - g, prime) == 0
    end
    println("Polynomial ModP Euclidean - PASSED")
end


#################################
##### OTHER GENERAL TESTING #####
#################################

function CRT_Prod_Test()
    N = 10
    for _ in 1:N
        p1 = rand(PolynomialSparseBInt)
        p2 = rand(PolynomialSparseBInt)

        prod = p1*p2
        comp_prod = CRT_prodBI(p1, p2)

        new2 = PolynomialSparseBInt([i for i in prod.terms])
        new1 = PolynomialSparseBInt([i for i in comp_prod.terms])

        test1 = degrees(new1)
        test2 = degrees(new2)
        test3 = coeffs(new1)
        test4 = coeffs(new2)

        @assert (test1[1] == test2[1] && test3[1] == test4[1])
    end
    println("CRT Test - Passed")

end


function PolyModP_generality_test()
    N = 100
    prime = 5
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        p2 = rand(PolynomialSparse)

        p1mod = PolynomialModP(p1, prime)
        p2mod = PolynomialModP(p2, prime)

        test1 = PolynomialSparse([i for i in mod(p1+p2, 5)])
        test2 = (p1mod + p2mod).p

        @assert degrees(test1)[1] == degrees(test2)[1]
        
    end

    println("PolynomialModP Generality - Passed")
end

function PolyModP_generality_test_prime()
    N = 1000
    prime = 5
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        p2 = rand(PolynomialSparse)
        
        p1mod = PolynomialModP(p1, 101)
        p2mod = PolynomialModP(p2, 101)
        
        test1 = divide(p1, p2)(101)
        _,r = test1
        test11 = PolynomialSparse([i for i in r])
        _,r2 = divide(p1mod, p2mod)
        
        @assert r2 == test11
    end
    println("PolynomialModP prime testing - PASSED")
end


function sqr_exp_test()
    N = 100
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        n = rand(1:8)

        exp1 = p1^n
        exp2 = repsqr_exp(p1, n)

        @assert coeffs(exp1)[1] == coeffs(exp1)[1]
    end
    println("Squared method test - PASSED")
end

sqr_exp_test()

println("Product type testing:")
prod_test_poly_sparse()
prod_test_poly_sparseBI()
prod_derivative_test_poly_modp()

println()
println("Derivative type testing:")
prod_derivative_test_poly_sparse()
prod_derivative_test_poly_sparseBI()
prod_derivative_test_poly_modp()

println()
println("Division type testing:")
division_test_poly()
division_test_poly_modp()
division_test_polyBI()

println()
println("Euclidean type testing:")
ext_euclid_test_poly()
#ext_euclid_test_polyBI()
ext_euclid_test_poly_modp()

println()
println("PolynomialModP Generality & Prime Tests:")
PolyModP_generality_test()
PolyModP_generality_test()

println()
println("CRT product output test:")
CRT_Prod_Test()
println()
println("Squared Power method test:")
sqr_exp_test()

