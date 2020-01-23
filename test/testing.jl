using GAPGroups
using Test

#to be modified once defined the new type
GG=GAPGroups.GAPGroup
GGE=GAPGroups.GAPGroupElem

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
