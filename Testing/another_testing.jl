function extended_euclid_alg(a::PolynomialSparse, b::PolynomialSparse, prime::Int)
    old_r, r = mod(a, prime), mod(b, prime)
    old_s, s = one(PolynomialSparse), zero(PolynomialSparse)
    old_t, t = zero(PolynomialSparse), one(PolynomialSparse)

    while !iszero(mod(r,prime))
        q = first(divide(old_r, r)(prime))
        old_r, r = r, mod(old_r - q*r, prime)
        old_s, s = s, mod(old_s - q*s, prime)
        old_t, t = t, mod(old_t - q*t, prime)
    end
    g, s, t = old_r, old_s, old_t

    @assert mod(s*a + t*b - g, prime) == 0
    return g, s, t  
end

a = rand(PolynomialSparse)
b = rand(PolynomialSparse)
prime = 101

old_r, r = mod(a, prime), mod(b, prime)
old_s, s = one(PolynomialSparse), zero(PolynomialSparse)
old_t, t = zero(PolynomialSparse), one(PolynomialSparse)

#while !iszero(mod(r,prime))
q = first(divide(old_r, r)(prime))
old_r, r = r, mod(old_r - q*r, prime)
old_s, s = s, mod(old_s - q*s, prime)
old_t, t = t, mod(old_t - q*t, prime)
#end
g, s, t = old_r, old_s, old_t


x = x_poly()
c = 93*x^4 + 85*x^2 + 88*x^1
d = 63*x^4 + 45*x^2 + 43*x^1
prime = 101

old_u, u = mod(c, prime), mod(d, prime)
old_v, v = one(Polynomial), zero(Polynomial)
old_w, w = zero(Polynomial), one(Polynomial)

#while !iszero(mod(r,prime))
z = first(divide(old_u, u)(prime))
old_u, u = u, mod(old_u - z*u, prime)
old_v, v = v, mod(old_v - z*v, prime)
old_w, w = w, mod(old_w - z*w, prime)
#end
g, s, t = old_r, old_s, old_t



q,z = Term(0,0), Term(0,0)
p = 101

x = x_poly()
num = 73*x^4 + 61*x^2 + 97*x^1
den = 58*x^3 + 92*x
f, g = mod(num,p), mod(den,p)
prev_degree = degree(f)
while degree(f) ≥ degree(g) 
    h = Polynomial( (leading(f) ÷ leading(g))(p) )  #syzergy 
    f = mod((f - h*g), p)
    q = mod((q + h), p)  
    prev_degree == degree(f) && break
    prev_degree = degree(f)
end
mod((num  - (q*g + f)),p)

#num1 = rand(PolynomialSparse)
#den1 = rand(PolynomialSparse)
j, k = mod(num1,p), mod(den1,p)
prev_degree = degree(j)
while degree(j) ≥ degree(k) 
    l = PolynomialSparse( (leading(j) ÷ leading(k))(p) )  #syzergy 
    j = mod((j - l*k), p)
    z = mod((z + l), p)  
    prev_degree == degree(f) && break
    prev_degree = degree(j)
end



function divide(num::PolynomialSparse, den::PolynomialSparse)
    function division_function(p::Int)
        f, g = mod(num,p), mod(den,p)
        degree(f) < degree(num) && return nothing 
        iszero(g) && throw(DivideError())
        q = Term(0,0)
        prev_degree = degree(f)
        while degree(f) ≥ degree(g) 
            h = PolynomialSparse( (leading(f) ÷ leading(g))(p) )  #syzergy 
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



function extended_euclid_alg(a::PolynomialSparse, b::PolynomialSparse, prime::Int)
    old_r, r = mod(a, prime), mod(b, prime)
    old_s, s = one(PolynomialSparse), zero(PolynomialSparse)
    old_t, t = zero(PolynomialSparse), one(PolynomialSparse)

    while !iszero(mod(r,prime))
        q = first(divide(old_r, r)(prime))
        old_r, r = r, mod(old_r - q*r, prime)
        old_s, s = s, mod(old_s - q*s, prime)
        old_t, t = t, mod(old_t - q*t, prime)
    end
    g, s, t = old_r, old_s, old_t

    @assert mod(s*a + t*b - g, prime) == 0
    return g, s, t  
end


a = rand(PolynomialSparse)
b = rand(PolynomialSparse)

x = x_poly()
aa = 54*x^5 + 99*x^4 + 6*x^3 + 74*x^1
bb = 72*x^4 + 79*x^3 + 26*x^2
extended_euclid_alg(aa, bb, 101)
extended_euclid_alg(a, b, 101)
prime = 101
old_r, r = mod(b, prime), mod(a, prime)
q = (divide(old_r, r)(prime))

divide(b, a)(101)


old_r, r = mod(a, prime), mod(b, prime)
old_s, s = one(PolynomialSparse), zero(PolynomialSparse)
old_t, t = zero(PolynomialSparse), one(PolynomialSparse)

while !iszero(mod(r,prime))
    q = first(divide(old_r, r)(prime))
    old_r, r = r, mod(old_r - q*r, prime)
    old_s, s = s, mod(old_s - q*s, prime)
    old_t, t = t, mod(old_t - q*t, prime)
end
g, s, t = old_r, old_s, old_t

@assert mod(s*a + t*b - g, prime) == 0
return g, s, t  
end