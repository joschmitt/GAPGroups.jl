

@testset "Cosets" begin
  
  G = dihedral_group(8)
  H, mH = center(G)
  
  @test index(G, H) == 4
  
  C = right_coset(H, G[1])
  @test order(C) == length(elements(C))
  
  @test length(right_cosets(G, H)) == index(G, H)
  
  @test length(right_transversal(G, H)) == index(G, H)

end

@testset "Conjugacy classes" begin
  
  G = symmetric_group(5)
  
  cc = conjugacy_class(G, G[1])
  cc1 = conjugacy_class(G, G[1]^G[2])
  @test order(cc) == order(cc1)
  @test cc == cc1
  

end
