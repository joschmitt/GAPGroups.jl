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
 #  @test typeof(G)==GG
end

@testset "Elements of Sym($i)" for i in 4:9
   G=symmetric_group(i)
   x=rand(G)
   y=rand(G)
   ox=order(x)
   oy=order(y)
   
#   @test x isa GGE
   @test ox isa Int64
   @test (x*y).X == x.X*y.X
   @test x^2 == x*x
   @test x*identity(G)==x
   @test x*x^-1 == identity(G)
   @test x^(ox-1)==x^-1
   @test x^y == y^-1*x*y
   @test x>y || x==y || x<y
end

end
