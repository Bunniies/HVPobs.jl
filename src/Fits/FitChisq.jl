function gen_uncorrelated_chisq(f::Function, x::Array{<:Real}, err::Vector{Float64}) #constrained
    chisq(par, dat) = sum((dat .- f(x,par)).^2 ./err.^2)
    return chisq
end

function gen_correlated_chisq(f::Function, x::Array{<:Real}, invCov::Matrix{Float64})
    chisq(par, dat) = begin
        r = dat .- f(x, par)
        v = invCov * r
        return sum(r .* v)
    end
    # (dat .- f(x,par))' * invCov * (dat .- f(x,par))
    return chisq
end