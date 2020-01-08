#testing GAPGroups.jl

include("GAPGroups.jl")

n=10
G=symmetric_group(n)

order(G) == factorial(n)

x=rand(G)
y=rand(G)

ox=order(x)
oy=order(y)

(x^y).X == (y^(oy-1)*x*y).X
