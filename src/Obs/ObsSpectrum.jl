@doc raw"""
    meff(obs::Vector{uwreal}, plat::Vector{Int64}; pl::Bool=false, data::Bool=false, wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)
    meff(corr::Corr, plat::Vector{Int64}; pl::Bool=true, data::Bool=false, wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)

Computes effective mass for a given correlator corr at a given plateau `plat`.
Correlator can be passed as an `Corr` struct or `Vector{uwreal}`.
    
The flags `pl` and `data` allow to show the plots and return data as an extra result.

```@example
data = read_mesons(path, "G5", "G5")
corr_pp = corr_obs.(data)
m = meff(corr_pp[1], [50, 60], pl=false)
```
"""
function meff(obs::Vector{uwreal}, plat::Vector{Int64}; pl::Bool=false, data::Bool=false, wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)

    tvals = length(obs)
    m = 0.5 .* log.((obs[2:tvals-2] ./ obs[3:tvals-1]) .^2)
    # m = log.((obs[2:tvals-2] ./ obs[3:tvals-1]) )

    m_av = plat_av(m, plat, wpm=wpm)

    if pl
        isnothing(wpm) ? uwerr(m_av) : uwerr(m_av, wpm)

        errorbar(collect(1:length(m)), value.(m), err.(m), fmt="d", mfc="none", color="black", capsize=2)
        fill_between(plat[1]:plat[2], value(m_av)-err(m_av), value(m_av)+err(m_av), alpha=0.6, color="royalblue")
        ylabel(L"$m_\mathrm{eff}$")
        xlabel(L"$x_0$")
        ylim(value(m_av)-80*err(m_av), value(m_av)+80*err(m_av))


        display(gcf())
        close("all")
    end

    !data ? (return m_av) : (return m_av, m)
end
meff(corr::Corr, plat::Vector{Int64}; pl::Bool=false, data::Bool=false, wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing) = meff(corr.obs, plat; pl=pl, data=data, wpm=wpm)


@doc raw"""
    mpcac(a0p::Vector{uwreal}, pp::Vector{uwreal}, plat::Vector{Int64}; ca::Float64=0.0, pl::Bool=true, data::Bool=false, wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)

    mpcac(a0p::Corr, pp::Corr, plat::Vector{Int64}; ca::Float64=0.0, pl::Bool=true, data::Bool=false, wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)

Computes the bare PCAC mass for a given correlator `a0p` and `pp` at a given plateau `plat`.
Correlator can be passed as an `Corr` struct or `Vector{uwreal}`.
    
The flags `pl` and `data` allow to show the plots and return data as an extra result. The `ca` variable allows to compute `mpcac` using
the improved axial current.

```@example
data_pp = read_mesons(path, "G5", "G5")
data_a0p = read_mesons(path, "G5", "G0G5")
corr_pp = corr_obs.(data_pp)
corr_a0p = corr_obs.(data_a0p)
m12 = mpcac(corr_a0p, corr_pp, [50, 60], pl=false)

p0 = 9.2056
p1 = -13.9847
g2 = 1.73410
ca = -0.006033 * g2 *( 1 + exp(p0 + p1/g2))

m12 = mpcac(corr_a0p, corr_pp, [50, 60], pl=false, ca=ca)
```
"""
function mpcac(a0p::Vector{uwreal}, pp::Vector{uwreal}, plat::Vector{Int64}; ca::Float64=0.0, pl::Bool=true, data::Bool=false, 
    kappa::Union{Vector{Float64}, Nothing}=nothing, mu::Union{Vector{Float64}, Nothing}=nothing, 
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)     
    
    corr_a0p = -a0p[1:end] 
    corr_pp = pp[1:end]

    if ca != 0.0
        der_a0p = (corr_a0p[3:end] .- corr_a0p[1:end-2]) / 2 
        der2_pp = (corr_pp[1:end-4] - 2*corr_pp[3:end-2] + corr_pp[5:end])/4
        der_a0p = der_a0p[2:end-1] + ca * der2_pp
    else
        der_a0p = (corr_a0p[4:end-1] .- corr_a0p[2:end-3]) / 2 
    end

    #der_a0p = (corr_a0p[3:end] .- corr_a0p[1:end-2]) / 2 
    #if ca != 0.0
    #    der2_pp = corr_pp[1:end-2] + corr_pp[3:end] - 2 * corr_pp[2:end-1]
    #    der_a0p = der_a0p + ca * der2_pp
    #end

    aux = der_a0p ./ ( 2. .* corr_pp[3:end-2])
    mass = plat_av(aux, plat, wpm=wpm)                                                                            
    if pl
        isnothing(wpm) ? uwerr(mass) : uwerr(mass, wpm)                       
        x = 1:length(aux)
        y = value.(aux)
        dy = err.(aux)
        v = value(mass)
        e = err(mass)
        
        figure()
        fill_between(plat[1]:plat[2], v-e, v+e, color="green", alpha=0.75)
        errorbar(x, y, dy, fmt="x", color="black")
        ylabel(L"$m_\mathrm{PCAC}$")
        xlabel(L"$x_0$")
        _, max_idx = findmax(value.(pp))
        ylim(v-20*e, v+20*e)
        #xlim(left=max_idx)
    
        if !isnothing(kappa)
            title(string(L"$\kappa_1 = $", kappa[1], L" $\kappa_2 = $", kappa[2]))
        end
        if !isnothing(mu)
            title(string(L"$\mu_1 = $", mu[1], L" $\mu_2 = $", mu[2]))
        end
        display(gcf())
    end                                                  
    if !data                   
        return mass                                        
    else               
        return (mass, aux)  
    end                    
end   

function mpcac(a0p::Corr, pp::Corr, plat::Vector{Int64}; ca::Float64=0.0, pl::Bool=true, data::Bool=false,
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)
    
    if a0p.kappa == pp.kappa || a0p.mu == pp.mu
        if a0p.mu == [0.0, 0.0]
            return mpcac(a0p.obs, pp.obs, plat, ca=ca, pl=pl, data=data, kappa=a0p.kappa, mu=nothing, wpm=wpm)
        else
            return mpcac(a0p.obs, pp.obs, plat, ca=ca, pl=pl, data=data, mu=a0p.mu, kappa=nothing, wpm=wpm)
        end
    else
        error("mu or kappa values does not match")
    end
end