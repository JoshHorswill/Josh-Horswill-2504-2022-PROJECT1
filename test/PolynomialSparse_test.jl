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
    println("prod_test_poly_sparse - PASSED")
end

prod_test_poly_sparse()

"""
Test derivative of polynomials (as well as product).
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
    println("prod_derivative_test_poly - PASSED")
end

prod_derivative_test_poly_sparse()


"""
Test division of polynomials modulo p.
"""
function division_test_poly(;prime::Int = 101, N::Int = 10^4, seed::Int = 0)
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
    println("division_test_poly - PASSED")
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
    println("ext_euclid_test_poly - PASSED")
end
ext_euclid_test_poly()

"""
Test product of polynomialsparse types
"""
function prod_test_poly_sparseBI(;N::Int = 10^3, N_prods::Int = 20, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparseBInt)
        p2 = rand(PolynomialSparseBInt)
        prod = p1*p2
        @assert p1.terms[1] * p2.terms[1] == prod.terms[2]
    end
    println("prod_test_poly_sparse - PASSED")
end
prod_test_poly_sparse()

"""
Test derivative of polynomials (as well as product).
"""
function prod_derivative_test_poly_sparseBI(;N::Int = 10^2,  seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparseBInt)
        p2 = rand(PolynomialSparseBInt)
        p1d = derivative(p1)
        p2d = derivative(p2)
        @assert (p1d*p2 + p2d*p1).terms[2] == (derivative(p1*p2)).terms[2]
    end
    println("prod_derivative_test_poly - PASSED")
end

prod_derivative_test_poly_sparse()

"""
Test division of polynomials modulo p.
"""
function division_test_polyBI(;prime::Int = big"101", N::Int = 10^4, seed::Int = 0)
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
        @assert iszero( mod(q*p2+r - p_prod, prime) )
    end
    println("division_test_poly - PASSED")
end
division_test_poly()

function ext_euclid_test_polyBI(;prime::Int=big"101", N::Int = 10^3, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparseBInt)
        p2 = rand(PolynomialSparseBInt)
        g, s, t = extended_euclid_alg(p1, p2, prime)
        @assert mod(s*p1 + t*p2 - g, prime) == 0
    end
    println("ext_euclid_test_poly - PASSED")
end
ext_euclid_test_poly()