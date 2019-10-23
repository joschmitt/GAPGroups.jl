function gap_header(type::UInt8, flags::UInt8, size::Int64)
    @assert size < 2^48
    ans = UInt64(0)
    ans += type
    ans = ans << 8
    ans += flags
    ans = ans << 48
    ans += size
    return ans
end

gap_header(size::Int64) = gap_header(0x08, 0x00, size)

struct GAPPerm
    mem::Vector{Int64}

    function GAPPerm(d::Vector{Int})
        n = length(d)
        mem = Vector{Int}(undef, n+2)
        mem[1] = gap_header(n)
        mem[2] = C_NULL # void pointer
        mem[3:end] .= d
        return new(mem)
    end

    function GAPPerm(n::Int)
        mem = Vector{Int64}(undef, n+2)
        mem[1] = gap_header(n)
        mem[2] = C_NULL # void pointer
        mem[3:end] .= 1:n
        return new(mem)
    end
end

function Base.getindex(p::GAPPerm, n::Integer)
    @boundscheck n > 0
    return (n > degree(p) ? n : p.mem[n+2])
end

function Base.setindex!(p::GAPPerm, v::Integer, n::Integer)
    @boundscheck n > 0
    @boundscheck n <= degree(p)
    return p.mem[n+2] = v
end

gap_header(p::GAPPerm) = UInt(p.mem[1])
gap_type(p::GAPPerm) = UInt8(gap_header(p) >> 56)
degree(p::GAPPerm) = Int(gap_header(p) << 16 >> 16)


### Generic stuff below

Base.iterate(p::GAPPerm, s=1) = (s > degree(p) ? nothing : (p[s], s+1))
Base.eltype(::Type{GAPPerm}) = Int64
Base.length(p::GAPPerm) = degree(p)
Base.size(p::GAPPerm) = (degree(p),)

Base.similar(p::GAPPerm) = GAPPerm(degree(p))

Base.firstindex(p::GAPPerm) = 1
Base.lastindex(p::GAPPerm) = degree(p)

function Base.:(==)(p::GAPPerm, q::GAPPerm)
    last_idx = max(degree(p), degree(q))
    for i in 1:last_idx
        p[i] == q[i] || return false
    end
    return true
end

Base.hash(p::GAPPerm, h::UInt) = hash(GAPPerm, hash(view(p, 1, degree(p)), h))

function mul!(out::GAPPerm, p::GAPPerm, q::GAPPerm)
    @boundscheck degree(out) >= max(degree(p), degree(q))
    out = (out === p || out === q ? similar(out) : out)

    @inbounds for i in 1:degree(out)
      out[i] = q[p[i]]
   end
   return out
end



function Generic.inv!(out::GAPPerm, p::GAPPerm)
    @boundscheck degree(p) == degree(out)
    out = (out === p ? similar(out) : out)

    @inbounds for i in 1:degree(p)
      out[p[i]] = i
   end
   return out
end

function Base.:*(p::GAPPerm, q::GAPPerm)
    out = (degree(p) > degree(q) ? similar(p) : similar(q))
    @inbounds out = mul!(out, p, q)
    return out
end

Base.inv(g::GAPPerm) = @inbounds inv!(similar(g), g)


using Random, AbstractAlgebra, Test

@testset "GAPPerms" begin

    let n = 51
        for _ in 1:100

            p = randperm(n)
            q = randperm(n)
            P = GAPPerm(p)
            Q = GAPPerm(q)

            @test gap_header(P) == gap_header(Q)
            @test degree(P) == degree(Q) == n
            @test gap_type(P) == 0x08

            @test collect(P) isa Vector{Int64}
            @test length(collect(P)) == degree(P)

            @test P*Q isa GAPPerm
            R = P*Q
            @test degree(R) == n
            @test gap_header(R) == gap_header(P)

            aaP = Perm(p)
            aaQ = Perm(q)

            @test aaP*aaQ == R

            @test inv(P)*P == Perm(n)
        end
    end
end
