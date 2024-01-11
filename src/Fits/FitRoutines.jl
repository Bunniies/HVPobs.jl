function fit_routine(model::Function, xdata::Array{<:Real}, ydata::Array{uwreal}, param::Int64=3; info::Bool=false, pval::Bool=false, wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)

    isnothing(wpm) ? uwerr.(ydata) : [uwerr(yaux, wpm) for yaux in ydata]
    
    yval = value.(ydata)
    yer = err.(ydata)
    
    @info("Uncorrelated fit")
    chisq = gen_uncorrelated_chisq(model, xdata, yer)
    fit = curve_fit(model, xdata, yval, 1.0 ./ yer.^2, fill(0.5, param))
    (upar, chi_exp) = isnothing(wpm) ? fit_error(chisq, coef(fit), ydata) : fit_error(chisq, coef(fit), ydata, wpm)
    chi2_fit_res = sum(fit.resid.^2 )
    
    # compute and print single point contribution to chi2
    # println("\n")
    # for i in 1:length(fit.resid)
        # println((fit.resid[i])^2)
    # end
    # println("\n")

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

    println("χ2/χ2exp: ", chi2_fit_res, " / ", chi_exp, " (dof: ", length(yval) - param,")")
    
    if !pval
        return FitRes(length(yval) - param, upar, chi2_fit_res, chi_exp, nothing)
    else
        pvalue = get_pvalue(chisq, chi2_fit_res, value.(upar), ydata, wpm=wpm, W= 1 ./ err.(ydata).^2, nmc=10000)
        println("p-value: ", pvalue)
        return FitRes(length(yval) - param, upar, chi2_fit_res, chi_exp, pvalue)        
    end

end

function fit_routine(model::Vector{Function}, xdata::Vector{Array{Float64, N}} where N, ydata::Vector{Array{uwreal, N}} where N, param::Int64; wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing,
    correlated_fit::Bool=false, pval::Bool=false)

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
        @info("Uncorrelated fit ")
        (upar, chi2_exp) = isnothing(wpm) ? fit_error(chisq, sol.minimizer, dat) : fit_error(chisq, sol.minimizer, dat, wpm)
    else
        @info("Correlated fit ")

        (upar, chi2_exp) = isnothing(wpm) ? fit_error(chisq, sol.minimizer, dat, W=inv_cov_tot) : fit_error(chisq, sol.minimizer, dat, wpm, W=inv_cov_tot)
    end

    println("χ2/χ2exp: ", min_fun(sol.minimizer), " / ", chi2_exp, " (dof: ", length(dat) - param,")")
    
    if !pval
        return FitRes(length(dat) - param, upar, min_fun(sol.minimizer), chi2_exp, nothing)
    else
        pvalue = get_pvalue(chisq, chi2_fit_res, value.(upar), dat, wpm=wpm, W= 1 ./ err.(dat).^2, nmc=10000)
        println("p-value: ", pvalue)
        return FitRes(length(yval) - data, upar, min_fun(sol.minimizer), chi2_exp, pvalue)        
    end
end