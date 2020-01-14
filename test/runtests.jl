using GAPGroups
using Test

#to be modified once defined the new type
GG=GAPGroups.GAPGroup
GGE=GAPGroups.GAPGroupElem

@testset "GAPGroups.jl" begin

n=10
G=symmetric_group(n)
@testset "The group Sym(n)" begin
   @test order(G) isa Int64
   @test order(G) == factorial(n)
   @test typeof(G)==GG
end

@testset "Elements of Sym($i)" for i in 4:9
   G=symmetric_group(i)
   x=rand(G)
   y=rand(G)
   z=perm([j for j=1:i])
   w=perm([1,2],[j for j=3:i])
   ox=order(x)
   oy=order(y)
   oz=order(z)

   @test x isa GGE
   @test ox isa Int64
   @test (x*y).X == x.X*y.X
   @test x^2 == x*x
   @test x*one(G)==x
   @test oz==i
   @test x*x^-1 == one(G)
   @test x^(ox-1)==x^-1
   @test x^w == w^-1*x*w
   @test z^y == y^-1*z*y
   @test z(1)==2
   @test (w^-1)(i)==i-1
   @test x(y(1))==(y*x)(1)
   @test x>y || x==y || x<y
   @test isequal(x,y) || isless(x,y) || isless(y,x)
   @test (isless(x,y) && x<y) || (isequal(x,y) && x==y) || (isless(y,x) && x>y)
end

end
