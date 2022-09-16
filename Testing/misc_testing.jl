function divide(num::Polynomial, den::Polynomial)
    function division_function(p::Int)
        f, g = mod(num,p), mod(den,p)
        degree(f) < degree(num) && return nothing 
        iszero(g) && throw(DivideError())
        q = Polynomial()
        prev_degree = degree(f)
        while degree(f) ≥ degree(g) 
            h = Polynomial( (leading(f) ÷ leading(g))(p) )  #syzergy 
            f = mod((f - h*g), p)
            q = mod((q + h), p)  
            println(f, q)
            prev_degree == degree(f) && break
            prev_degree = degree(f)
        end
        #@assert iszero( mod((num  - (q*g + f)),p))
        return q, f
    end
    return division_function
end

p1 = Polynomial([Term(1,4), Term(1,5)])
p2 = Polynomial(Term(1,2))

divide(p1, p2)(101)


function divide(num::PolynomialSparse, den::PolynomialSparse)
    function division_function(p::Int)
        f, g = mod(num,p), mod(den,p)
        degree(f) < degree(num) && return nothing 
        iszero(g) && throw(DivideError())
        q = PolynomialSparse()
        prev_degree = degree(f)
        while degree(f) ≥ degree(g) 
            h = PolynomialSparse( (leading(f) ÷ leading(g))(p) )  #syzergy 
            f = mod((f - h*g), p)
            q = mod((q + h), p)  
            println(f, q)
            prev_degree == degree(f) && break
            prev_degree = degree(f)
        end
        #@assert iszero( mod((num  - (q*g + f)),p))
        return q, f
    end
    return division_function
end


function -(p::PolynomialSparse, t::Term)::PolynomialSparse
    pp = deepcopy(p)
    cnt = 0
    for terms in pp
        cnt += 1
        if terms.degree == t.degree
            pp.terms[cnt] += (t)
        else
            pp.terms[cnt] = terms
        end
    end

    if t.degree ∉ degrees(p)
        push!(pp, (t))
    end
    return pp
end


function -(p1::PolynomialSparse, p2::PolynomialSparse)::PolynomialSparse
    p = deepcopy(p1)
    for t in p2
        p = p + t
    end
    return p
end


p3 = PolynomialSparse([Term(1,1), Term(10, 3), Term(4, 8), Term(2,2)])
p3^2


p4 = PolynomialSparse([Term(1,1), Term(10, 2)])
t = Term(2, 9)
p = deepcopy(p3)

if t.degree ∉ degrees(p)
    push!(p, t)
end
p



p = PolynomialSparse([Term(1,2), Term(2,3), Term(5,9)])
t = Term(5,2)

function +(p::PolynomialSparse, t::Term)
    cnt = 0
    pp = deepcopy(p)

    if t.degree ∉ degrees(pp)
        push!(pp, t)
    end

    for i in pp
        cnt += 1
        if i.degree == t.degree
            pp.terms[cnt] += t
        else
            pp.terms[cnt] = i
        end
    end
    return pp
end

p + t

p = PolynomialSparse([Term(3,2), Term(2,3), Term(5,9)])
t = Term(3,2)
function -(p::PolynomialSparse, t::Term)
    cnt = 0
    pp = deepcopy(p)

    if t.degree ∉ degrees(pp)
        push!(pp, (-t))
    end

    for i in pp
        cnt += 1
        if i.degree == t.degree
            pp.terms[cnt] += (-t)
        else
            pp.terms[cnt] = i
        end
    end
    return pp
end
x = p - t


function +(p1::PolynomialSparse, p2::PolynomialSparse)::PolynomialSparse
    p = deepcopy(p1)
    for t in p2
        p = p + t
    end
    return p
end

x = PolynomialSparse([Term(3,2), Term(2,3), Term(5,9), Term(5,10)])
y = PolynomialSparse([Term(3,2), Term(2,3)])
x + y


function -(p1::PolynomialSparse, p2::PolynomialSparse)::PolynomialSparse
    p = deepcopy(p1)
    for t in p2
        p = p - t
    end
    return p
end

x = PolynomialSparse([Term(3,2), Term(2,3), Term(5,9), Term(5,10)])
y = PolynomialSparse([Term(3,2), Term(2,3), Term(5,9), Term(5,10)])
z = x - y
typeof(z)

p = 101
num = Polynomial(Term(1,4))
den = Polynomial(Term(1,2))
f, g = mod(num,p), mod(den,p)

num = PolynomialSparse(Term(1,4))
den = PolynomialSparse(Term(1,2))
j, k = mod(num,p), mod(den,p)

q1 = Polynomial()
prev_degree = degree(f)
h1 = Polynomial( (leading(f) ÷ leading(g))(p) )  #syzergy 
f = mod((f - h1*g), p)
q1 = mod((q1 + h1), p)

q2 = PolynomialSparse()
prev_degree = degree(j)
h2 = PolynomialSparse( (leading(j) ÷ leading(k))(p) )  #syzergy 
j = mod((j - h2*k), p)
q2 = mod((q2 + h2), p)


h1.terms[3] == h2.terms[1]

f
h1*g
h1
g


f - h1*g

j
h2
g
h2*k
