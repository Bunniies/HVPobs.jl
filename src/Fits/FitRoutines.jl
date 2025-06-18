@doc raw"""
    fit_routine(model::Function, xdata::Array{<:Real}, ydata::Array{uwreal}, param::Int64=3; wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)
    fit_routine(model::Vector{Function}, xdata::Vector{Array{Float64, N}} where N, ydata::Vector{Array{uwreal, N}} where N, param::Int64; wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing, correlated_fit::Bool=false, pval::Bool=false)


Given a `model` function with a number `param` of parameters,  an array of `xdata` and an array of `ydata` of type uwreal,
this function performs the fit with the given `model` and returns a `FitRes` object.

If an array of `model` is passed instead, together with a vector of array of `xdata` and `ydata`, the function performs a global fit by applying each 
model element-wise to each set of data.

This function automatically performs error propagation from `ydata` to the fit parameters returned in `FitRes`.

The flag `pval` allows to compute the fit pvalues, if set to true.

The flag `wpm` is an optional array of Float64 of lenght 4 used by `ADerrors` package for error analysis. The first three paramenters specify the criteria to determine
the summation windows:

- `vp[1]`: The autocorrelation function is summed up to ``t = round(vp[1])``.

- `vp[2]`: The sumation window is determined using U. Wolff poposal with ``S_\tau = wpm[2]``

- `vp[3]`: The autocorrelation function ``\Gamma(t)`` is summed up a point where its error ``\delta\Gamma(t)`` is a factor `vp[3]` times larger than the signal.

An additional fourth parameter `vp[4]`, tells ADerrors to add a tail to the error with ``\tau_{exp} = wpm[4]``.
Negative values of `wpm[1:4]` are ignored and only one component of `wpm[1:3]` needs to be positive.
If the flag `covar`is set to true, `fit_routine` takes into account covariances between x and y for each data point.
```@example
@. model(x,p) = p[1] + p[2] * exp(-(p[3]-p[1])*x)
@. model2(x,p) = p[1] + p[2] * x[:, 1] + (p[3] + p[4] * x[:, 1]) * x[:, 2]

# single fit
fit_routine(model, xdata, ydata, param=3)
# global fit
fit_routine([model, model2], [xdata, xdata2], [ydata, ydata2], param=3)
```
"""
function fit_routine(model::Function, xdata::Array{<:Real}, ydata::Array{uwreal}, param::Int64=3; pval::Bool=false, wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing, info::Bool=false, lineprint::Bool=true, fitRes::Bool=false, ensList::Union{Vector{String}, Nothing}=nothing)

    isnothing(wpm) ? uwerr.(ydata) : [uwerr(yaux, wpm) for yaux in ydata]
    
    yval = value.(ydata)
    yer = err.(ydata)
    
    lineprint ? @info("Uncorrelated fit") : nothing
    chisq = gen_uncorrelated_chisq(model, xdata, yer)
    fit = curve_fit(model, xdata, yval, 1.0 ./ yer.^2, fill(0.5, param))
    # println("chi2 ", sum(fit.resid.^2))
    # println("coef ", coef(fit))
    (upar, chi_exp) = isnothing(wpm) ? fit_error(chisq, coef(fit), ydata) : fit_error(chisq, coef(fit), ydata, wpm)
    chi2_fit_res = sum(fit.resid.^2 )
    
    # compute and print single point contribution to chi2
    if !isnothing(ensList) && lineprint
        println("Single ensemble contribution to chi2")
        for i in 1:length(fit.resid)
            println("  $(ensList[i]) => $((fit.resid[i])^2)")
        end
        println("------------------------------------")
    end

    if !isapprox(chi2_fit_res, value(chisq(coef(fit), ydata)), atol=1e-2)
       @warn "Chi2 from the two determination not within 0.01 tolerance "
        println("chi2 from fit residual:    ", chi2_fit_res)
        println("chi2 from custom function: ", chisq(coef(fit), ydata))
    end

    if info
        for i in eachindex(upar)
            isnothing(wpm) ? uwerr(upar[i]) : uwerr(upar[i], wpm)
            print("\n Fit parameter: ", i, ": ")
            details(upar[i])
        end
    end

    lineprint ? println("χ2/χ2exp: ", chi2_fit_res, " / ", chi_exp, " (dof: ", length(yval) - param,")") : nothing
    
    if !pval
        fitRes ? (return [FitRes(length(yval) - param, upar, chi2_fit_res, chi_exp, nothing),fit.resid]) : (return FitRes(length(yval) - param, upar, chi2_fit_res, chi_exp, nothing))
    else
        pvalue = get_pvalue(chisq, chi2_fit_res, value.(upar), ydata, wpm=wpm, W= 1 ./ err.(ydata).^2, nmc=10000)
        lineprint ? println("p-value: ", pvalue) : nothing
        fitRes ? ([FitRes(length(yval) - param, upar, chi2_fit_res, chi_exp, pvalue),fit.resid]) : (return FitRes(length(yval) - param, upar, chi2_fit_res, chi_exp, pvalue))
    end
end

function fit_routine(model::Vector{Function}, xdata::Vector{Array{Float64, N}} where N, ydata::Vector{Array{uwreal, N}} where N, param::Int64; wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing, correlated_fit::Bool=false, pval::Bool=false, lineprint::Bool=true)

    if !(length(model) == length(xdata) == length(ydata))
        error("Dimension mismatch")
    end
    N = length(model)

    dat = ydata[1]
    idx = Vector{Vector{Int64}}(undef, N)
    e = Vector{Vector{Float64}}(undef, N)
    j = 1
    for i = 1:N
        if isnothing(wpm)
            uwerr.(ydata[i])
        else
            for yaux in ydata[i]
                [uwerr(yaux[k], wpm) for k in eachindex(yaux)]
            end
        end
        
        e[i] = err.(ydata[i])
        if i > 1 
            dat = vcat(dat, ydata[i])
        end
        stp = j + length(ydata[i]) - 1
        idx[i] = collect(j:stp)
        j = stp + 1
    end
    if !correlated_fit
        chisq = (par, dat) -> sum([sum((dat[idx[i]] .- model[i](xdata[i], par)).^2 ./e[i].^2) for i=1:N])
    else 
        error("correlated fit not implemented")
    end

    min_fun(t) = chisq(t, value.(dat))
    p = fill(0.5, param)
    #println(chisq(p,dat))
    sol = optimize(min_fun, p, LBFGS()) 
    #println(chisq(sol.minimizer,dat))
    
    if !correlated_fit
        lineprint ? @info("Uncorrelated fit ") : nothing
        (upar, chi2_exp) = isnothing(wpm) ? fit_error(chisq, sol.minimizer, dat) : fit_error(chisq, sol.minimizer, dat, wpm)
    else
        lineprint ? @info("Correlated fit ") : nothing

        (upar, chi2_exp) = isnothing(wpm) ? fit_error(chisq, sol.minimizer, dat, W=inv_cov_tot) : fit_error(chisq, sol.minimizer, dat, wpm, W=inv_cov_tot)
    end

    lineprint ? println("χ2/χ2exp: ", min_fun(sol.minimizer), " / ", chi2_exp, " (dof: ", length(dat) - param,")") : nothing
    
    if !pval
        return FitRes(length(dat) - param, upar, min_fun(sol.minimizer), chi2_exp, nothing)
    else
        pvalue = get_pvalue(chisq, chi2_fit_res, value.(upar), dat, wpm=wpm, W= 1 ./ err.(dat).^2, nmc=10000)
        lineprint ? println("p-value: ", pvalue) : nothing
        return FitRes(length(yval) - data, upar, min_fun(sol.minimizer), chi2_exp, pvalue)        
    end
end


function lin_fit(x::Vector{<:Real}, v::Vector{Float64}, e::Vector{Float64})
    sig2 = e .* e
    S = sum(1 ./ sig2)
    Sx = sum(x ./ sig2)
    Sy = sum(v ./ sig2)
    Sxy = sum(v .* x ./ sig2)
    Sxx = sum(x .* x ./sig2)
    delta = S * Sxx - Sx*Sx
    par = [Sxx*Sy-Sx*Sxy, S*Sxy-Sx*Sy] ./delta
    #C = [[Sxx/delta, -Sx/delta], [-Sx/delta,  S/delta]]
    return par
end

@doc raw"""
    lin_fit(x::Vector{<:Real}, y::Vector{uwreal})

Computes a linear fit of uwreal data points y. This method return uwreal fit parameters and chisqexpected.

```@example
fitp, csqexp = lin_fit(phi2, m2)
m2_phys = fitp[1] + fitp[2] * phi2_phys
```
"""
function lin_fit(x::Vector{<:Real}, y::Vector{uwreal}; wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)
    isnothing(wpm) ? uwerr.(y) : [uwerr(yaux, wpm) for yaux in y]
    par = lin_fit(x, value.(y), err.(y))
    chisq(p, d) = sum((d .- p[1] .- p[2].*x).^2 ./ err.(y) .^2)
    (fitp, csqexp) = fit_error(chisq, par, y)
    for i in eachindex(fitp) 
        isnothing(wpm) ? uwerr(fitp[i]) : uwerr(fitp[i], wpm)
        print("\n Fit parameter: ", i, ": ")
        details(fitp[i])
    end
    println("Chisq / chiexp: ", chisq(par, y), " / ", csqexp, " (dof: ", length(x)-length(par),")")
    return (fitp, csqexp)
end

@doc raw"""
    x_lin_fit(par::Vector{uwreal}, y::Union{uwreal, Float64})

Computes the results of a linear interpolation/extrapolation in the x axis
"""
x_lin_fit(par::Vector{uwreal}, y::Union{uwreal, Float64}) = (y - par[1]) / par[2]
@doc raw"""
    y_lin_fit(par::Vector{uwreal}, y::Union{uwreal, Float64})

Computes the results of a linear interpolation/extrapolation in the y axis
"""
y_lin_fit(par::Vector{uwreal}, x::Union{uwreal, Float64}) = par[1] + par[2] * x






# FIT ROUTINES FROM A RAMOS

function fit_defs_yerr(f, x, dy)
    
    function lmfit(prm, y)
        
	nof = round(Int64, length(y))
	res = Vector{eltype(prm)}(undef,nof)
	# for i in 1:nof
	    # res[i] = (y[i] - f(x[i], prm)) / dy[i]
	# end
    res =   (y .- f(x, prm)) ./ dy
	return res
    end
    chisq(prm,data) = sum(lmfit(prm, data) .^ 2)
    
    return lmfit, chisq
end

function fit_defs_yerr_corr(f, x, dy, Winv)

    u = LinearAlgebra.cholesky(LinearAlgebra.Symmetric(Winv)).U

    function lmfit(prm, y)
        
	nof  = round(Int64, length(y))
	res  = Vector{eltype(prm)}(undef,nof)
	res2 = Vector{eltype(prm)}(undef,nof)
    res =   (y .- f(x, prm)) ./ dy

	# for i in 1:nof
	    # res[i] = (y[i] - f(x[i], prm))
	# end

        for k in 1:length(res)
            res2[k] = u[k,1]*res[1]
            for i in 2:length(res)
                res2[k] = res2[k] + u[k,i]*res[i]
            end
        end

	return res2
    end
    chisq(prm,data) = sum(lmfit(prm, data) .^ 2)
    
    return lmfit, chisq
end

function fit_data_yerr(f, xv, yv, nparam; correlated=false)

    ADerrors.uwerr.(yv)
    if correlated
        cvi = LinearAlgebra.inv(ADerrors.cov(yv))
        lm, csq = fit_defs_yerr_corr(f, xv, ADerrors.err.(yv), cvi)

        prm0 = fill(0.5,nparam) #zeros(nparam)
        fit  = LeastSquaresOptim.optimize(xx -> lm(xx, ADerrors.value.(yv)), prm0,
		                          LeastSquaresOptim.LevenbergMarquardt(), autodiff = :forward)
    
        fitp, csqexp = ADerrors.fit_error(csq, fit.minimizer, yv, W=cvi)
    else
        lm, csq = fit_defs_yerr(f, xv, ADerrors.err.(yv))
        prm0 = zeros(nparam)
        fit  = LeastSquaresOptim.optimize(xx -> lm(xx, ADerrors.value.(yv)), prm0,
		                          LeastSquaresOptim.LevenbergMarquardt(), autodiff = :forward)

        fitp, csqexp = ADerrors.fit_error(csq, fit.minimizer, yv)
    end
    println("χ2/χ2exp: ", fit.ssr, " / ", csqexp, " (dof: ", length(yv) - nparam,")")

    return FitRes(length(yv)-nparam, fitp, fit.ssr, csqexp, nothing)
    # return fitp, csqexp, fit.ssr    
end

function fit_defs_xyerr(f, dxy, Nalpha)
    
    function lmfit(prm, xy)
        Ndata = div(length(xy), Nalpha+1)
        Npar = length(prm) - Ndata*Nalpha
        # println("Nalpha: ", Nalpha)
        # println("Ndata:  ", Ndata)
        # println("Npar:   ", Npar)
        p = prm[1:Npar]
        res = Vector{Vector{eltype(prm)}}(undef, Ndata)
        for k = 1:Ndata
            inverr = [1 / dxy[k + (i-1)*Ndata] for i=1:Nalpha+1]
            xx = [prm[Npar + k + (i-1)*Ndata] for i=1:Nalpha]

            delta = [xy[k + (i-1)*Ndata] - xx[i] for i=1:Nalpha]
            yy = f(xx', p)
            push!(delta, xy[k+Nalpha*Ndata] - yy[1])
            res[k] = delta .* inverr
        end

        return vcat(res...)
	   
    end

    chisq(prm,data) = sum(lmfit(prm, data) .^ 2)
    
    return lmfit, chisq
end

function fit_data_xyerr(f, xv, yv, nparam)
    Nalpha = size(xv,2) # number of x-variables
    Ndata = size(yv, 1) # number of datapoints
    ADerrors.uwerr.(xv)
    ADerrors.uwerr.(yv)

    dxyv = [vcat(value.(xv)...);value.(yv)]
    dxye = [vcat(err.(xv)...);err.(yv)]
    dxy = [vcat(xv...);yv]

    lm, csq = fit_defs_xyerr(f, dxye, Nalpha)
    
    prm0 = zeros(nparam+length(xv))
    fit  = LeastSquaresOptim.optimize(xx -> lm(xx, dxyv), prm0,
		    LeastSquaresOptim.LevenbergMarquardt(), autodiff = :forward)
    
    fitp, csqexp = ADerrors.fit_error(csq, fit.minimizer, dxy)
    
    println("χ2/χ2exp: ", fit.ssr, " / ", csqexp, " (dof: ", length(yv) - nparam,")")
    return FitRes(length(yv) - nparam, fitp, fit.ssr, csqexp, nothing )
    # return fitp, csqexp, fit.ssr
end

function fit_data(f, xv, yv, nparam; correlated=false)
    if (isa(xv, Array{ADerrors.uwreal}))
        return fit_data_xyerr(f, xv, yv, nparam)
    else
        return fit_data_yerr(f, xv, yv, nparam, correlated=correlated)
    end
end
