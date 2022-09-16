function factor(f::PolynomialSparse, prime::Int)::Vector{Tuple{PolynomialSparse,Int}}
    #Cantor Zassenhaus factorization

    f_modp = mod(f, prime)
    degree(f_modp) ≤ 1 && return [(f_modp,1)]

    # make f primitive
    ff = prim_part(f_modp)(prime)      
    # @show "after prim:", ff

     # make f square-free
    squares_poly = gcd(f, derivative(ff), prime) 
    ff = (ff ÷ squares_poly)(prime) 
    # @show "after square free:", ff

    # make f monic
    old_coeff = leading(ff).coeff
    ff = (ff ÷ old_coeff)(prime)        
    # @show "after monic:", ff

    dds = dd_factor(ff, prime)

    ret_val = Tuple{PolynomialSparse,Int}[]

    for (k,dd) in enumerate(dds)
        sp = dd_split(dd, k, prime)
        sp = map((p)->(p ÷ leading(p).coeff)(prime),sp) #makes the polynomials inside the list sp, monic
        for mp in sp
            push!(ret_val, (mp, multiplicity(f_modp,mp,prime)) )
        end
    end

    #Append the leading coefficient as well
    push!(ret_val, (leading(f_modp).coeff* one(PolynomialSparse), 1) )

    return ret_val
end


x = PolynomialSparse([Term(1,2), Term(3,4)])
factor(x, 3)

x = Polynomial([Term(1,2), Term(3,4)])
factor(x,3)