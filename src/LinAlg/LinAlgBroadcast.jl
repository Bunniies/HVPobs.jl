

#Base.:*(a::uwreal, b::AbstractArray) = fill(a, size(b)) .* b

function Base.iterate(uw::uwreal)
    return iterate(uw, 1)
end

function Base.iterate(uw::uwreal, state::Int64)
    if state > length(uw)
        return nothing
    else
        return uw[state], state +1
    end
end

function Base.iterate(uw::Vector{uwreal})
    return iterate(uw, 1)
end

function Base.iterate(uw::Vector{uwreal}, state::Int64)
    if state > length(uw)
        return nothing
    else
        return uw[state], state +1
    end
end

function Base.length(uw::uwreal)
    return length(value(uw))
end

function Base.getindex(uw::uwreal, ii::Int64)
    idx = getindex(value(uw), ii)
    return  uw
end