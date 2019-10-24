using GAPGroups
using AbstractAlgebra
using Random
using Test

@testset "GAPGroups.jl" begin

    import GAPGroups.GAPPerm

    Base.:(==)(aap::Generic.Perm, gapp::GAPPerm) = aap.d == collect(gapp)

    @testset "GAPPerms" begin
        n = 129
        p = randperm(n)
        q = randperm(n)
        P = GAPPerm(p)
        Q = GAPPerm(q)
        aaP = Perm(p)
        aaQ = Perm(q)

        @testset "gap_header" begin
            import GAPGroups: gap_header, gap_flags, gap_type

            (deg, fl, tp) = (2^32+2, UInt8(9), UInt8(1))
            gh = gap_header(deg, fl, tp)
            tmp = similar(Q);
            tmp.data[1:2] .= reinterpret(UInt32, [gh])

            @test degree(tmp) == deg
            @test gap_flags(tmp) == fl
            @test gap_type(tmp) == tp

            @test gap_header(P) == gap_header(Q)
            @test degree(P) == degree(Q) == n
            @test gap_type(P) == 0x08
            @test gap_flags(P) == 0x00
        end

        @testset "GAPperms operations" begin
            @test collect(P) isa Vector{Int64}
            @test length(collect(P)) == degree(P)

            @test P*Q isa GAPPerm
            R = P*Q
            @test degree(R) == n
            @test gap_header(R) == gap_header(P)

            @test aaP*aaQ == R
            @test inv(P)*P == GAPPerm(n)

            S = GAPPerm(randperm(100))
            @test P*S isa GAPPerm
            @test S*P isa GAPPerm

            @test degree(P*S) == degree(S*P) == 129 # degree(P)

            perms = [randperm(129) for _ in 1:100]

            @test all(all(p .== GAPPerm(p)) for p in perms)

            @test all(Perm(p)*Perm(q) == GAPPerm(p)*GAPPerm(q)
                for p in perms, q in perms)

            perms = [randperm(rand(1:1000)) for _ in 1:100]

            @test all(degree(GAPPerm(p)*GAPPerm(q)) == max(length(p), length(q))
                for p in perms, q in perms)
        end
    end
end
