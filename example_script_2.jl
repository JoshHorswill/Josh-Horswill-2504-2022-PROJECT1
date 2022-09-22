using Pkg
Pkg.activate(".")

#include("poly_factorization_project.jl")

# Construct Polynomial
x = x_poly()

# Declare two functions
x = x_poly()
p1 = 5x^3 + 8x^2 - 9x
p2 = 1x^2 - 2x + 4x + 3 

# Print functions
println("Using x_poly, we construst two polynomial functions:")
println("   1) ", p1)
println("   2) ", p2)

println() #Multiplication
println("These functions are operable:")
println("Multiplication p1 × p2 =")
println(p1*p2)

println() #Division
div = (p1÷p2)(101);
println("Division (w.r.t prime 101) p1 ÷ p2 = $div")

y = PolynomialSparse([Term(1,1), Term(-4, 3), Term(4, 0)])
y1 = PolynomialSparse([Term(1,3), Term(2, 2)])

println() #Addition
println("Addition p1 + p2 =")
println(p1+p2)

println() #Subtraction
println("Subtraction p1 - p2 =")
println(p1-p2)

# Derivative example
println()
println("The polynomial functions are also able to be derived via the derivative implementation:")
println("Polynomial Function: ", p1)
println("→ Derivative: ", derivative(p1))
println()
println("Where the derivative function can be considered as an operator:")
println("Operation considered: p1(x) × p2'(x)")
println("Result = ", p1*derivative(p2))

# Modulo example
prime = 19
println()
println("We can find the mod result of a function to a given prime:")
println("Prime = ", prime)
println("Function = ", p1)
println("function mod prime = ",  mod(p1,prime))

# GCD
println()
println("As an extension to modulo, we can find the greatest common denominator of two polynomial functions:")
println("Function 1: ", p1, " & Function 2: ", p2)
println("Greatest Common Denominator: ", gcd(p1, p2, 19))


