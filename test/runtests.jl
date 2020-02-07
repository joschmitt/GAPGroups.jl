using GAPGroups
using Test

#to be modified once defined the new type
GG=GAPGroups.GAPGroup
GGE=GAPGroups.GAPGroupElem


function test_properties(n::Int64)

   G= @inferred GG symmetric_group(n)
   @testset "The group Sym(n)" begin
      @test deg(G) isa Integer
      @test deg(G) == n
      @test hasorder(G) isa Bool
      @test hasorder(G)
      @test order(G) isa Int64
      @test order(Int32, G) isa Int32
      @test order(G) == factorial(n)
      @test order(BigInt, G) == factorial(n)
      @test typeof(G)==GG
      @test hasgens(G) isa Bool
      @test hasgens(G)
      @test gens(G) isa Vector{GGE}
   end
end

function explicit_example(n::Int64)

   if n>9
   @testset "Explicit example" begin
      x=cperm(1:5,6:8,9:n)
      A=vcat([i for i in 10:n],[9])
      y=perm(vcat([2,3,4,5,1,7,8,6],A))

      @test parent(x)==symmetric_group(n)
      @test x==y
      @test order(x) == lcm(15,n-8)
      @test x(3)==4
      @test x(8)==6
      @test x(n)==9
      H=symmetric_group(n+3)
      @test parent(x)!=H
      @test parent(H(x))==H
      x=H(x)
      @test parent(x)==H
   end
   end
end

function test_operations(L::Union{Array{Int64,1},UnitRange{Int64}})
   @testset "Elements of Sym($i)" for i in L
      if i>1
      G=symmetric_group(i)
      x=@inferred GGE rand(G)
      y=rand(G)
      #z=cperm(1:i)
      z=perm(vcat([j for j in 2:i],[1]))
      if minimum(L)>3
         w=cperm([1,2],[j for j=3:i])
      else
         w=cperm([1,2])
      end
      ox= @inferred Int64 order(x)
      oy=order(y)
      oz=order(z)

      @test x isa GGE
      @test ox isa Int64
      @test sign(z)==(-1)^(i-1)
      @test sign(x*y)==sign(x)*sign(y)
      @test parent(x)==G
      @test parent(z)==G
      @test x/y == x*y^-1
      @test (x*y).X == x.X*y.X
      @test x^2 == x*x
      @test x*one(G)==x
      @test oz==i
      @test x*x^-1 == one(G)
      @test x^(ox-1)==x^-1
      @test x^w == w^-1*x*w
      @test z^y == y^-1*z*y
      @test conj(x,y) == y^-1*x*y
      @test comm(x,y) == x^-1*conj(x,y)
      @test comm(x,y) == div_left(x^y,x)
      @test div_right(x,y) == x/y
      @test z(1)==2
      @test (z^-1)(i)==i-1
      @test x(y(1))==(y*x)(1)
      @test x>y || x==y || x<y
      @test isequal(x,y) || isless(x,y) || isless(y,x)
      @test (isless(x,y) && x<y) || (isequal(x,y) && x==y) || (isless(y,x) && x>y)
      end
   end
end

