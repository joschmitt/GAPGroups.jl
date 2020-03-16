

@testset "Cosets" begin
  
   G = dihedral_group(8)
   H, mH = center(G)
  
   @test index(G, H) == 4
  
   C = right_coset(H, G[1])
   @test order(C) == length(elements(C))
  
   @test length(right_cosets(G, H)) == index(G, H)
  
   @test length(right_transversal(G, H)) == index(G, H)

   @testset "set comparation for cosets in PermGroup" begin
      G=symmetric_group(5)
      x = G(cperm([1,2,3]))
      y = G(cperm([1,4,5]))
      z = G(cperm([1,2],[3,4]))
      H = sub(G,[y])[1]
      K = sub(G,[z])[1]

      @test_throws ArgumentError rc = right_coset(H,K(z))
      @test_throws ArgumentError rc = left_coset(H,K(z))
      rc = right_coset(H,x)
      lc = left_coset(H,x)
      dc = double_coset(H,x,K)
      @test acting_domain(rc) == H
      @test representative(rc) == x
      @test representative(lc) == x
      @test Set(elements(rc)) == Set([h*x for h in H])
      @test Set(elements(lc)) == Set([x*h for h in H])
      @test Set(elements(dc)) == Set([h*x*k for h in H for k in K])
   end

end

@testset "Conjugacy classes" begin
  
  G = symmetric_group(5)
  
  cc = conjugacy_class(G, G[1])
  cc1 = conjugacy_class(G, G[1]^G[2])
  @test order(cc) == order(cc1)
  @test cc == cc1
  

end
