@doc raw"""
    bayesian_av(fun::Function, y::Array{uwreal}, tmin_array::Array{Int64}, tmax_array::Array{Int64}, k::Int64, pl::Bool, data::Bool; wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)

    bayesian_av(fun1::Function, fun2::Function, y::Array{uwreal}, tmin_array::Array{Int64}, tmax_array::Array{Int64}, k1::Int64, k2::Int64, pl::Bool, data::Bool; wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)

    bayesian_av(fun::Array{Function}, y::Array{uwreal}, tmin_array::Array{Int64}, tmax_array::Array{Int64}, k::Array{Int64}, pl::Bool, data::Bool; wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)

Computes bayesian average of data. For a given fit function, it explores choices of fit intervals, assigning each of them a weight. The function saves the `first` fit parameter of your function, and then it does the weighted average of it and assigns a systematic. See https://arxiv.org/abs/2008.01069 

The function takes as input the fit intervals to explore. 

`tmin_array` is an array of integers with the lower bounds on the fit intervals to explore, ***ordered from lower to higher***.

`tmax_array` is an array of integers with the upper bounds on the fit intervals to explore, ***ordered from lower to higher***.

`k` is the number of parameters of the fit function to use.

You can also use as input two fit functions, and two values of `k`, one for each function. Then, for each fit interval choice, the function explores the two fit functions. This means that for each fit interval choice you get two results: one for the first fit funcction, and another for the second. You can also use a vector of functions and a vector of k (numer of parameters of each funtion) to apply the bayesian averaging method to multiple functions.

The method returns two objects: first, the weighted average as an uwreal object, with mean value and statistichal error. The second object returned is the systematic error coming from the fit interval variation. If `data` is `true`, then returns 4 objects: weighted average, systematic error, a vector with the results of the fit for each fit interval choice, and a vector with the weights associated to each fit.

```@example
@.fun(x,p) = p[1] * x ^0
k = 1
tmin_array = [10,11,12,13,14,15]
tmax_array = [80,81,82,83,84,85]
(average, systematics, data, weights) = bayesian_av(fun,x,tmin_array,tmax_array,k,pl=true,data=true)

@.fun1(x,p) = p[1] * x ^0
@.fun2(x,p) = p[1] + p[2] * exp( - p[3] * (x))
k1 = 1
k2 = 3
tmin_array = [10,11,12,13,14,15]
tmax_array = [80,81,82,83,84,85]
(average, systematics) = bayesian_av(fun1,fun2,x,tmin_array,tmax_array,k1,k2)

@.fun1(x,p) = p[1] * x ^0
@.fun2(x,p) = p[1] + p[2] * exp( - p[3] * (x))
k1 = 1
k2 = 3
tmin_array = [10,11,12,13,14,15]
tmax_array = [80,81,82,83,84,85]
(average, systematics) = bayesian_av([fun1,fun2],x,tmin_array,tmax_array,[k1,k2])
```
"""
function bayesian_av(fun::Function, y::Array{uwreal}, tmin_array::Array{Int64}, tmax_array::Array{Int64}, k::Int64; pl::Bool=false, data::Bool=false,
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing, path_plt::Union{String,Nothing}=nothing, 
    plt_title::Union{Nothing,String}=nothing, label::Union{Nothing, LaTeXString}=nothing, pbc::Bool=false)
    
    weight_model = Array{Float64,1}()
    AIC = Array{Float64,1}()
    chi2chi2exp = Array{Float64,1}()
    p1 = Array{uwreal,1}()
    mods = Array{String,1}()
    Pval = Array{Float64,1}()

    if tmax_array[end] > length(y)
        error("Error: upper bound for the fits is bigger than last data point")
    end

    total = length(y)
    isnothing(wpm) ? uwerr.(y) : [uwerr(y[i],wpm) for i in 1:length(y)]
    
    for INDEX in tmin_array ## vary tmin
        for j in tmax_array ## vary tmax
            if pbc==true
                #println("pbc is true, T= $(INDEX + j), y = $(length(y))")
                if INDEX + j != length(y) 
                    continue
                else
                   #println("if you see me: $(INDEX + j)==$(length(y))")
                end
            end
            #try
                x = [i for i in INDEX:1:j] .-1	
	            yy = y[INDEX:1:j]
	            Ncut = total - length(x)
                dy = err.(yy)
                W = 1 ./ dy .^2

                p00 = [0.5 for i in 1:1:k]
                chisq = gen_uncorrelated_chisq(fun, x, dy)

                fit = curve_fit(fun, x, value.(yy), W, p00)
                uwerr.(yy)
                isnothing(wpm) ? (up, chi_exp) = fit_error(chisq,coef(fit),yy) : (up,chi_exp) = fit_error(chisq,coef(fit),yy,wpm)
                isnothing(wpm) ? uwerr(up[1]) : uwerr(up[1],wpm)
                chi2 = sum(fit.resid.^2) 
                Q = get_pvalue(chisq, sum(fit.resid.^2), value.(up), yy, wpm=wpm, W=W, nmc=10000)
                
                push!(Pval, Q)

                # push!(AIC, chi2 + 2*k + 2*Ncut)   # AIC criteria
                push!(AIC, chi2 - 2*chi_exp)        # TIC criteria

                push!(chi2chi2exp, chi2 / chi_exp) 
                push!(p1, up[1])
                push!(mods,string("[", INDEX, ",", j, "]"))
            #catch 
            #   @warn string(":/ Negative window for error propagation at tmin = ", INDEX, ", tmax = ", j, "; skipping that point")
            #end
        end
    end
    
    # compute all weights
    offset = minimum(AIC)
    AIC = AIC .- offset
    weight_model = exp.(-0.5 .* AIC)
    weight_model = weight_model ./ sum(weight_model)
    
    # sort weights and discard 5% tail 
    idxW = sortperm(weight_model, rev=true)
    cumulative_w = cumsum(weight_model[idxW])
    # println("acceptance in bayeasian_av set to 0.99, should be 0.95")
    idxcumw = findfirst(x->x>=0.99, cumulative_w)
    # println(length(idxW) - idxcumw)
    idxW = sort(idxW[1:idxcumw])
    
    # take only the accepted 
    AIC = AIC[idxW]
    Pval = Pval[idxW]
    chi2chi2exp  = chi2chi2exp[idxW]
    mods = mods[idxW]
    p1 = p1[idxW]
    weight_model = weight_model[idxW]
    weight_model ./= sum(weight_model)


    
    p1_mean = sum(p1 .* weight_model)
    isnothing(wpm) ? uwerr(p1_mean) : uwerr(p1_mean, wpm) 
    # println(sum(p1 .^ 2 .* weight_model))
    # println((sum(p1 .* weight_model)) ^ 2)
    systematic_err =  0.0 #sqrt(sum(p1 .^ 2 .* weight_model) - (sum(p1 .* weight_model)) ^ 2) 
    if pl   
        fig = figure(figsize=(10,7.5))

        subplots_adjust(hspace=0.1) 
        subplot(411)    
        if !isnothing(plt_title)
            #title(plt_title)
        end
        ax1 = gca()                
        x = 1:length(p1)
        y = value.(p1)
        dy = err.(p1)
        v = value(p1_mean)
        e = err(p1_mean)

        fill_between(x, v-e, v+e, color="red", alpha=0.2)
        errorbar(x, y, dy, fmt="^", mfc="none", color="black", ms=10, capsize=2)
        setp(ax1.get_xticklabels(),visible=false) # Disable x tick labels
        isnothing(label) ?  ylabel(L"$p_1$") : ylabel(label)
    
        subplot(412)
        ax2=gca()
        bar(x, weight_model, alpha=0.4, color="royalblue", edgecolor="blue", linewidth=1.5)
        setp(ax2.get_xticklabels(),visible=false) # Disable x tick labels
        # errorbar(mods, weight_model, 0*dy, color="green")
        ylabel(L"$W$")

        subplot(413)
        ax3=gca()
        setp(ax3.get_xticklabels(),visible=false) # Disable x tick labels

        bar(x, Pval, alpha=0.4, color="royalblue", edgecolor="blue", linewidth=1.5 )
        ylabel(L"$p-value$")
        
        
        subplot(414)
        bar(x, chi2chi2exp, alpha=0.4, color="royalblue", edgecolor="blue", linewidth=1.5 )
        ylabel(L"$\chi^2/\chi^2_{\mathrm{exp}}$")
        xticks(x, mods, rotation=45)
        xlabel(L"$\mathrm{Models} \ [t_{\mathrm{min}}/a,t_{\mathrm{max}}/a]$")
        tight_layout()
        display(fig)
        if !isnothing(path_plt)
            #tt = plt_title *".pdf" 
            savefig(path_plt)
        end
        close("all")
    end

    if !data                   
        return (p1_mean, systematic_err)                                         
    else               
        FitStorage = Dict(
        "W"        => weight_model,
        "AIC"      => AIC .+ offset,
        "pval"     => Pval, 
        "chivsexp" => chi2chi2exp,
        "mods"     => mods
        )
        return (p1_mean, systematic_err, p1, FitStorage)    
    end 

end