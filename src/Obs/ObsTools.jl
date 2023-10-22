function plat_av(obs::Vector{uwreal}, plat::Vector{Int64}; wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)
    isnothing(wpm) ? uwerr.(obs) : [uwerr(obs[k], wpm) for k in eachindex(obs)]
    w = 1 ./ err.(obs[plat[1]:plat[2]]).^2
    av = sum(w .* obs[plat[1]:plat[2]]) / sum(w)
    return av
end