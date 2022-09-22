ps1 = PolynomialSparse([Term(4294967296,2)])
ps2 = PolynomialSparse([Term(4294967296, 2)])

ps2 + 10000000000000000*ps1

pb1 = PolynomialSparse([Term(4294967296,2)])
pb2 = PolynomialSparse([Term(4294967296, 2)])

pb2 + 10000000000000000*ps1

x = Term(922337203685477580, 1)

println("Range limits for Int64:")
println("    -> -9,223,372,036,854,775,808")
println("    ->  9,223,372,036,854,775,807")

t = Term(9223372036854775807, 1)
println()

println("Maximum Term coeff: $t")
println()
println("Overflow Example: Adding 1x^1 to the maximum representation:")
x = Term(1,1) + t
println("    -> $x")
println()
println("If we use BigInt, we can represent the same term:")
t1 = TermBI(big"9223372036854775807",1)
y = t1 + TermBI(big"1", 1)
println("     -> $y")
println()
println("We can go even further and represent a much larger range of Integers:")
t3 = TermBI(big"92233720392233720368547758076854775807", 1)
println("    -> $t3")
println()

t4 = TermBI(big"200", 1)*t3
println("And perform arithmatic:")
println("    -> 200 * Term = $t4")

println()
println("################ TIMING ################")
println("Comparision of overflow times")

pls = PolynomialSparse(Term(9223372036854775807, 2))
plt = PolynomialSparse(Term(25, 2))

plsBI = PolynomialSparseBInt(TermBI(big"9223372036854775807", 2))
pltBI = PolynomialSparseBInt(TermBI(big"25", 2))

println()
println("Addition of $pls and $plt")
println("PolynomialSparse:")
@time pls + plt
println("PolynomialSparseBInt:")
@time plsBI + pltBI 

println()
println("Product of $pls and $plt")
println("PolynomialSparse:")
@time pls * plt
println("PolynomialSparseBInt:")
@time plsBI * pltBI 

println()
println("Comparision of compiling 10x^1 + 20x^2 + 100x^10:")

println("Polynomial")
@time PolynomialDense([Term(10,1), Term(20, 2), Term(100, 10)])
println("PolynomialSparse")
@time PolynomialSparse([Term(10,1), Term(20, 2), Term(100, 10)])
println("PolynomialSparseBInt")
@time PolynomialSparseBInt([TermBI(big"10",1), TermBI(big"20", 2), TermBI(big"100", 10)])

println()
println("Operation Timing: product of two random polynomials")
println("Polynomial")
@time rand(PolynomialDense)*rand(PolynomialDense)
println("PolynomialSparse")
@time rand(PolynomialSparse)*rand(PolynomialSparse)
println("PolynomialSparseBInt")
@time rand(PolynomialSparseBInt)*rand(PolynomialSparseBInt)

