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