@doc raw"""
    pvalue(chisq::Function,
    chi2::Float64,
    xp::Vector{Float64}, 
    data::Vector{uwreal},
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing} = Dict{Int64,Vector{Float64}}();
    W::Union{Vector{Float64},Array{Float64,2}} = Vector{Float64}(),
    nmc::Int64 = 5000)

Computes the p-value of a previously done fit, using as input the `\chi^2` observed from the fit, the fit parameters and the fitted data. 
The p-value for a given `\chi^2` is the probability of, given the data you have, finding such a `\chi^2` or worse from a fit, and still
have the data well described by the fit function. `nmc` is the number of MC samples used to estimate the p-value integral, default is 5000.
By now it only works with a vector for weights (containing the diagonal of W)

```@example
function fit_defs(f::Function,x,W)
    chisq(p,d)=sum((d-f(x,p)).^2 .*W)
    return chisq
end

@.fun(x,p) = p[1] + p[2] * x
chisq = fit_defs(fun, x, 1.0 ./ err.(y) .^ 2)
fit = curve_fit(fun, x, value.(y), 1.0 ./ err.(y) .^ 2, [0.5, 0.5])
(up, chiexp) = fit_error(chisq, coef(fit), y)

wpm = Dict{Int64,Vector{Float64}}()
wpm[1] = [-1.0,-1.0,4-0,-1.0]
Q = pvalue(chisq, chi2, value.(up), y, wpm; W = 1.0 ./ err.(y) .^ 2, nmc=10000)
#Q = pvalue(chisq, chi2, value.(up), y; W = 1.0 ./ err.(y) .^ 2, nmc=10000)
#Q = pvalue(chisq, chi2, value.(up), y)
```
"""
function get_pvalue(chisq::Function,
    chi2::Float64,
    xp::Vector{Float64}, 
    data::Vector{uwreal};
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing} =  Dict{Int64,Vector{Float64}}(),
    W::Union{Vector{Float64},Array{Float64,2}} = Vector{Float64}(),
    nmc::Int64 = 5000)

    n = length(xp)   # Number of fit parameters
    m = length(data) # Number of data

    xav = zeros(Float64, n+m)
    for i in 1:n
        xav[i] = xp[i]
    end
    for i in n+1:n+m
        xav[i] = data[i-n].mean
    end
    ccsq(x::Vector) = chisq(view(x, 1:n), view(x, n+1:n+m)) 
    if (n+m < 4)
        cfg = ForwardDiff.HessianConfig(ccsq, xav, ADerrors.Chunk{1}());
    else
        cfg = ForwardDiff.HessianConfig(ccsq, xav, ADerrors.Chunk{4}());
    end
        
    hess = Array{Float64}(undef, n+m, n+m)
    ForwardDiff.hessian!(hess, ccsq, xav, cfg)
        
    cse = 0.0
    Q = dQ = 0.0
    if (m-n > 0)
        if (length(W) == 0)
            Ww = zeros(Float64, m)
            for i in 1:m
                if (data[i].err == 0.0)
                    #isnothing(wpm) ? wuerr(data[i]) : uwerr(data[i], wpm)
                    uwerr(data[i], wpm)
                    if (data[i].err == 0.0)
                        error("Zero error in fit data")
                    end
                end
                Ww[i] = 1.0 / data[i].err^2
            end
        else
            Ww = W
        end
        #cse = chiexp(hess, data, Ww, wpm)

        m = length(data)
        n = size(hess, 1) - m

        hm = view(hess, 1:n, n+1:n+m)
        sm = Array{Float64, 2}(undef, n, m)

        for i in 1:n, j in 1:m
            sm[i,j] = hm[i,j] / sqrt.(Ww[j])
        end
        maux = sm * sm'
        hi   = LinearAlgebra.pinv(maux)
        Px   = -hm' * hi * hm

        for i in 1:m
            Px[i,i] = Ww[i] + Px[i,i]
        end
        C = cov(data) 
        
        nu = sqrt(C) * Px * sqrt(C)
        
        N = length(nu[1,:])
        z = randn(N, nmc)

        eig = abs.(eigvals(nu))
        eps = 1e-14 * maximum(eig)
        eig = eig .* (eig .> eps)

        aux = eig' * (z .^ 2)
        Q = 1.0 - mean(aux .< chi2)

        x = chi2 .- eig[2:end]' * (z[2:end,:].^2)
        x = x / eig[1]
        #dQ = juobs.mean((x .> 0) .* exp.(-x * 0.5) * 0.5 ./ sqrt.(abs.(x)))
        #dQ = err(cse)/value(cse) * dQ
    end

    return Q #uwreal([Q,dQ],"")
end