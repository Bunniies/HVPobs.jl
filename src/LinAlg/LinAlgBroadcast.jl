

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


function Base.:*(x::uwreal, y::Array{Any})
    N = size(y, 1)
    return fill(x, N) .* y
end  
Base.:*(x::Array{Any}, y::uwreal) = Base.:*(y,x)

function Base.:*(x::uwreal, y::Array{Float64})
    N = size(y, 1)
    return fill(x, N) .* y
end
Base.:*(x::Array{Float64}, y::uwreal) = Base.:*(y,x)

function Base.:*(x::uwreal, y::Array{uwreal})
    N = size(y, 1)
    return fill(x, N) .* y
end
Base.:*(x::Array{uwreal}, y::uwreal) = Base.:*(y,x)

function Base.:+(x::uwreal, y::Vector{uwreal})
    N = size(y, 1)
    return fill(x, N) .+ y
end

# overloading uwreal * Rational
function Base.:*(x::uwreal, y::Rational{T}) where {T<:Real}
    num = y.num * x
    return Rationaluw(num, uwreal(Float64(y.den)))
end
Base.:*(y::Rational{T}, x::uwreal) where {T<:Real}  = Base.:*(x,y)

# overloading Rationaluw * Rationaluw
function Base.:*(x::Rationaluw, y::Rationaluw)
    num = x.num * y.num
    den = x.den * y.den
    return Rationaluw(num, den)
end
