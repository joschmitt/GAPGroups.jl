using GAPGroups
using Test


function test_properties(n::Int64)

   G= @inferred PermGroup symmetric_group(n)
   @testset "The group Sym($n)" begin
      @test degree(G) isa Integer
      @test degree(G) == n
      @test isfinite(G)
      @test order(G) isa Int64
      if n < 13 
        @test order(Int32, G) isa Int32 
      end
      @test order(G) == factorial(n)
      @test order(BigInt, G) == factorial(n)
      @test gens(G) isa Vector{PermGroupElem}
   end
end

function test_iterate(n::Int64)
  if n > 3 && n < 11          # otherwise the group is too big
    @testset "Iteration" begin
      G = symmetric_group(n)
      L = elements(G)
      @test L isa Vector{PermGroupElem}
      @test length(L) == factorial(G.deg)
      @test length(unique(L)) == factorial(G.deg)
      @test rand(G) isa PermGroupElem
      @test rand(G) in G
      A = PermGroup[]
      for x in G
        push!(A, x)
      end
      @test length(A) == factorial(G.deg)
      @test length(unique(A)) == factorial(G.deg)
      s = 0         # check if the number of (n-1)-cycles is correct
      for x in G 
        if order(x) == (n-1)
          s+=1
        end
      end
      @test s == factorial(n-2)*n
    end
  end
end

function explicit_example(n::Int64)

  if n > 9
    G=symmetric_group(n)
    @testset "Explicit example" begin
      x=cperm(1:5,6:8,9:n)
      A=vcat([i for i in 10:n],[9])
      A=vcat([2,3,4,5,1,7,8,6],A)
      y=perm(A)

      @test x==y
      @test A==listperm(y)
      @test x==perm(listperm(x))
      @test order(x) == lcm(15,n-8)
      @test x(3)==4
      @test x(8)==6
      @test x(n)==9
   end
   @testset "Change of parent" begin
      x=cperm(1:5,6:8,9:n)
      H=symmetric_group(n+3)
      K=symmetric_group(n-1)
      
      @test parent(x)==G
      x=cperm(H,1:5,6:8,9:n)
      @test parent(x)===H
      x=cperm(G,1:5,6:8,9:n)
      @test_throws ArgumentError x=cperm(K,1:5,6:8,9:n)
      @test parent(x)===G
      @test parent(H(x))==H
      @test parent(H(x))===H
      @test parent(x)!=H
      x=H(x)
      @test parent(x)===H
      @test_throws ArgumentError K(x)
   end
   end
end

function test_operations(L::Union{Array{Int64,1},UnitRange{Int64}})
   @testset "Elements of Sym($i)" for i in L
      if i>1
      G=symmetric_group(i)
      x=@inferred PermGroupElem rand(G)
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

      @test x isa PermGroupElem
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
      @test y==perm(listperm(y))
      end
   end
end

