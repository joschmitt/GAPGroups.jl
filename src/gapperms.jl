function gap_header(size::Int64, flags::UInt8=0x00, type::UInt8=0x08)
    @assert 0 < size < 2^48
    header = UInt64(0)
    header += size
    header = header << 16
    header += flags
    header = header << 8
    header += type
    return header
end

struct GAPPerm
    data::Vector{UInt32}

    function GAPPerm(d::AbstractVector{<:Integer})
        n = length(d)
        data = Vector{UInt32}(undef, n+4)
        data[1:2] = reinterpret(UInt32, [gap_header(n)])
        data[3:4] = reinterpret(UInt32, [C_NULL]) # void pointer
        data[5:end] .= convert.(UInt32, d)
        return new(data)
    end

    function GAPPerm(n::Int)
        data = Vector{UInt64}(undef, n+4)
        data[1:2] = reinterpret(UInt32, [gap_header(n)])
        data[3:4] = reinterpret(UInt32, [C_NULL]) # void pointer
        data[5:end] .= UInt32(1):UInt32(n)
        return new(data)
    end
end

function Base.getindex(p::GAPPerm, n::Integer)
    @boundscheck 0 < n
    return (n > degree(p) ? Int(n) : Int(p.data[n+4]))
end

function Base.setindex!(p::GAPPerm, v::Integer, n::Integer)
    @boundscheck 0 < n <= degree(p)
    return p.data[n+4] = v
end

gap_header(p::GAPPerm) = reinterpret(UInt64, p.data[1:2])[]
gap_flags(p::GAPPerm) = reinterpret(UInt8, p.data[1:1])[2]
gap_type(p::GAPPerm)  = reinterpret(UInt8, p.data[1:1])[1]
Generic.degree(p::GAPPerm) = Int(gap_header(p) >> 24)

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

function Generic.mul!(out::GAPPerm, p::GAPPerm, q::GAPPerm)
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
