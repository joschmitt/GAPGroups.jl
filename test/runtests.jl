using GAPGroups
using Test

#to be modified once defined the new type
GG=GAPGroups.GAPGroup
GGE=GAPGroups.GAPGroupElem

#=
@testset "GAPGroups.jl" begin

n=10
G= @inferred GG symmetric_group(n)
@testset "The group Sym(n)" begin
   @test order(G) isa Int64
   @test order(G) == factorial(n)
   @test typeof(G)==GG
end

@testset "Explicit example" begin
   x=perm([1,2,3,4,5],[6,7,8],[9,10,11,12])

   @test order(x) == 60
   @test x(3)==4
   @test x(8)==6
end

@testset "Elements of Sym($i)" for i in 4:9
   G=symmetric_group(i)
   x=@inferred GGE rand(G)
   y=rand(G)
   z=perm([j for j=1:i])
   w=perm([1,2],[j for j=3:i])
   ox= @inferred Int64 order(x)
   oy=order(y)
   oz=order(z)

   @test x isa GGE
   @test ox isa Int64
   @test sign(w)==(-1)^(i-2)
   @test sign(x*y)==sign(x)*sign(y)
   @test x/y == x*y^-1
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
=#

function test_GAPGroups_interface_conformance(g::GAPGroups.GAPGroupElem, h::GAPGroups.GAPGroupElem)
   #to be modified once defined the new type
   GG=GAPGroups.GAPGroup
   GGE=GAPGroups.GAPGroupElem

   @test parent(g) isa GG
   @test parent(h) isa GG
   G, H = parent(g), parent(h)

   @testset "Parent methods" begin
      @test elem_type(G) == typeof(g)
      @test one(G) isa typeof(g)
      @test one(G)==one(g)==one(h)

      @test hasorder(G) isa Bool
      @test hasgens(G) isa Bool
      @test ngens(G) isa Integer
      @test gens(G) isa Vector{typeof(g)}
      if hasorder(G)
         @test order(G) isa Integer
         @test order(G) > 0
      end
   end

   @testset "Comparison methods" begin
      @test (g==h) isa Bool
      @test isequal(g,h) isa Bool
      @test g == g
      @test isequal(h,h)
      @test (g<h) isa Bool
      @test isless(g,h) isa Bool
      @test g>h || g==h || g<h
      @test isequal(g,h) || isless(g,h) || isless(h,g)
   end

   @testset "Group operations" begin
      g1,h1 = deepcopy(g), deepcopy(h)

      @test inv(g) isa typeof(g)
      @test (g,h) == (g1,h1)
      @test g*h isa typeof(g)
      @test (g,h) == (g1,h1)
      @test g^2 == g*g
      @test (g,h) == (g1,h1)
      @test g^-3 == inv(g)*inv(g)*inv(g)
      @test (g,h) == (g1,h1)
      @test (g*h)^-1 == inv(h)*inv(g)
      @test (g,h) == (g1,h1)
      @test conj(g,h) == inv(h)*g*h
      @test (g,h) == (g1,h1)
      @test comm(g,h) == g^-1*h^-1*g*h
      @test (g,h) == (g1,h1)
      @test isone(g*inv(g)) && isone(inv(g)*g)
   end

   @testset "In-place operations" begin
      g1,h1 = deepcopy(g), deepcopy(h)
      out = rand(G)  #to be replaced by out=similar(g)

      @test isone(one!(g))
      g = deepcopy(g1)

      @testset "mul!" begin
         @test mul!(out,g,h) == g1*h1
         @test (g,h) == (g1,h1)

         @test mul!(out,g,h) == g1*h1
         @test (g,h) == (g1,h1)

         @test mul!(g,g,h) == g1*h1
         @test h==h1
         g = deepcopy(g1)

         @test mul!(h,g,h) == g1*h1
         @test g == g1
         h = deepcopy(h1)

         @test mul!(g,g,g) == g1*g1
         g = deepcopy(g1)
      end

      @testset "conj!" begin
         res = h1^-1*g1*h1
         @test conj!(out,g,h) == res
         @test (g,h) == (g1,h1)

         @test conj!(g,g,h) == res
         @test h == h1
         g = deepcopy(g1)

         @test conj!(h,g,h) == res
         @test g == g1
         h = deepcopy(h1)

         @test conj!(g,g,g) == g1
         g = deepcopy(g1)
      end

      @testset "comm!" begin
         res = g1^-1*h1^-1*g*h

         @test comm!(out,g,h) == res
         @test (g,h) == (g1,h1)

         @test comm!(g,g,h) == res
         @test h == h1
         g = deepcopy(g1)

         @test comm!(h,g,h) == res
         @test g == g1
         h = deepcopy(h1)
      end

      @testset "div_[left|right]!" begin
         res = g*h^-1

         @test div_right!(out,g,h) == res
         @test (g,h) == (g1,h1)

         @test div_right!(g,g,h) == res
         @test h == h1
         g = deepcopy(g1)

         @test div_right!(h,g,h) == res
         @test g == g1
         h = deepcopy(h1)

         @test div_right!(g,g,g) == one(g)
         g = deepcopy(g1)

         res = h^-1*g

         @test div_left!(out,g,h) == res
         @test (g,h) == (g1,h1)

         @test div_left!(g,g,h) == res
         @test h == h1
         g = deepcopy(g1)

         @test div_left!(h,g,h) == res
         @test g == g1
         h = deepcopy(h1)

         @test div_left!(g,g,g) == one(g)
         g = deepcopy(g1)
      end
   end
end
