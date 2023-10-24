function gen_uncorrelated_chisq(f::Function, x::Array{<:Real}, err::Vector{Float64}) #constrained
    chisq(par, dat) = sum((dat .- f(x,par)).^2 ./err.^2)
    return chisq
end