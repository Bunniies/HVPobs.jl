function meff(obs::Vector{uwreal}, plat::Vector{Int64}; pl::Bool=false, data::Bool=false, wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)

    tvals = length(obs)
    m = 0.5 .* log.((obs[2:tvals-2] ./ obs[3:tvals-1]) .^2)

    m_av = plat_av(m, plat)

    if pl
        isnothing(wpm) ? uwerr(m_av) : uwerr(m_av, wpm)

        errorbar(collect(1:length(m)), value.(m), err.(m), fmt="d", mfc="none", color="black", capsize=2)
        fill_between(plat[1]:plat[2], value(m_av)-err(m_av), value(m_av)+err(m_av), alpha=0.6, color="royalblue")
        ylabel(L"$m_\mathrm{eff}$")
        xlabel(L"$x_0$")

        display(gcf())
        close()
    end

    !data ? (return m_av) : (return m_av, m)
end