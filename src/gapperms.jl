struct GAPPerm{T<:Union{UInt16, UInt32}}
    data::Vector{T}

    function GAPPerm{T}(itr) where T<:Union{UInt16, UInt32}
        @boundscheck !(itr isa Integer) || throw(ArgumentError("Pass range 1:$itr to the GAPPerm{T} constructor."))
        n = length(itr)

        @assert 0 < n <= typemax(T)
        gaph_s = div(sizeof(UInt64), sizeof(T)) # 2 or 4
        ptr_s = div(sizeof(Ptr), sizeof(T)) # 2 or 4
        # size of GAPPerm data in bytes, including the pointer to inverse perm:
        size = T(n * sizeof(T) + sizeof(Ptr))

        data = Vector{T}(undef, n + gaph_s + ptr_s)
        data[1 : gaph_s] = reinterpret(T, [perm_header(size)])
        data[1+gaph_s : gaph_s+ptr_s] = reinterpret(T, [C_NULL]) # void pointer

        data[1+data_offset(GAPPerm{T}):end] .= itr

        return new(data)
    end
end

GAPPerm(itr) = GAPPerm{UInt32}(itr)
GAPPerm(n::Integer) = GAPPerm{UInt32}(1:n)

function gap_header(size::Integer, flags::UInt8, type::UInt8)
    @assert 0 <= size <= 2^48
    return (UInt64(size)<<16) | (UInt64(flags) << 8) | UInt64(type)
end

perm_header(size_inbytes::T) where T<:Union{UInt16, UInt32} =
    gap_header(size_inbytes, 0x00, (sizeof(T) == 4 ? 0x08 : 0x07))

data_offset(::Type{GAPPerm{T}}) where T = div(sizeof(UInt64) + sizeof(Ptr), sizeof(T))

Base.@propagate_inbounds function Base.getindex(p::T, n::Integer) where T<:GAPPerm
    @boundscheck 0 < n
    return (n > degree(p) ? Int(n) : Int(p.data[n+data_offset(T)]))
end

Base.@propagate_inbounds function Base.setindex!(p::T, v::Integer, n::Integer) where T<:GAPPerm
    @boundscheck 0 < n <= degree(p)
    return p.data[n+data_offset(T)] = v
end

function gap_header(p::GAPPerm{T}) where T
    gaph_s = div(sizeof(UInt64), sizeof(T))
    return reinterpret(UInt64, view(p.data, 1:gaph_s))[1]
end

gap_flags(p::GAPPerm) = reinterpret(UInt8, view(p.data, 1:1))[2]
gap_type(p::GAPPerm)  = reinterpret(UInt8, view(p.data, 1:1))[1]
gap_size(p::GAPPerm) = gap_header(p) >> 16

Generic.degree(p::GAPPerm{T}) where T = Int(div(gap_size(p) - sizeof(Ptr), sizeof(T)))

### Generic stuff below

Base.iterate(p::GAPPerm, s=1) = (s > degree(p) ? nothing : (p[s], s+1))
Base.eltype(::Type{<:GAPPerm}) = Int64
Base.length(p::GAPPerm) = degree(p)
Base.size(p::GAPPerm) = (degree(p),)

Base.similar(p::GAPPerm{T}) where T = GAPPerm{T}(1:degree(p))

Base.@propagate_inbounds function Base.:(==)(p::GAPPerm, q::GAPPerm)
    last_idx = max(degree(p), degree(q))
    @inbounds for i in 1:last_idx
        p[i] == q[i] || return false
    end
    return true
end

Base.hash(p::GAPPerm, h::UInt) = hash(GAPPerm, hash(view(p, 1, degree(p)), h))

Base.@propagate_inbounds function Generic.mul!(out::GAPPerm, p::GAPPerm, q::GAPPerm)
    @boundscheck degree(out) >= max(degree(p), degree(q))
    out = (out === p || out === q ? similar(out) : out)

    @inbounds for i in 1:degree(out)
      out[i] = q[p[i]]
   end
   return out
end

Base.@propagate_inbounds function Generic.inv!(out::GAPPerm, p::GAPPerm)
    @boundscheck degree(p) == degree(out)
    out = (out === p ? similar(out) : out)

    @inbounds for i in 1:degree(p)
      out[p[i]] = i
   end
   return out
end

Base.@propagate_inbounds function Base.:*(p::GAPPerm, q::GAPPerm)
    out = (degree(p) >= degree(q) ? similar(p) : similar(q))
    @inbounds out = mul!(out, p, q)
    return out
end

Base.@propagate_inbounds Base.inv(g::GAPPerm) = @inbounds inv!(similar(g), g)
