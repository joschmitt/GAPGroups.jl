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
        P = GAPPerm{UInt32}(p)
        Q = GAPPerm{UInt16}(q)
        aaP = Perm(p)
        aaQ = Perm(q)

        @testset "gap_header" begin
            import GAPGroups: gap_header, gap_flags, gap_type, gap_size

            @test gap_header(P) ==  0x00000000020c0008
            @test gap_size(P) == UInt(129*sizeof(UInt32) + sizeof(Ptr))
            @test degree(P) == 129
            @test gap_flags(P) == 0x00
            @test gap_type(P) == 0x08

            @test gap_header(Q) == 0x00000000010a0007
            @test gap_size(Q) == UInt(129*sizeof(UInt16) + sizeof(Ptr))
            @test degree(Q) == 129
            @test gap_flags(Q) == 0x00
            @test gap_type(Q) == 0x07

            for p in [P, Q]
                R = similar(p)
                @test degree(R) == degree(p)
                @test gap_flags(R) == gap_flags(p)
                @test gap_type(R) == gap_type(p)
                @test gap_header(R) == gap_header(p)
            end

            s = 2^32

            (deg, fl, tp) = (UInt(4*s+8), UInt8(9), UInt8(1))
            gh = gap_header(deg, fl, tp)
            tmp = similar(P);
            tmp.data[1:2] .= reinterpret(UInt32, [gh])

            @test degree(tmp) == s
            @test gap_flags(tmp) == fl
            @test gap_type(tmp) == tp
        end

        @testset "GAPperms operations" begin
            @test collect(P) isa Vector{Int}
            @test length(collect(P)) == degree(P)

            @test P*Q isa GAPPerm
            R = P*Q
            @test degree(R) == n
            @test gap_header(R) == gap_header(P)

            @test aaP*aaQ == R
            @test inv(P)*P == GAPPerm(n)

            S = GAPPerm{UInt32}(randperm(100))
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
