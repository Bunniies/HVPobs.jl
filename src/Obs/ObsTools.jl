function apply_rw(data::Array{Float64}, W::Matrix{Float64}, vcfg::Union{Nothing, Vector{Int64}}=nothing)
    nc =  isnothing(vcfg) ? collect(1:size(data, 1)) : vcfg
    W1 = W[1, nc]
    W2 = W[2, nc]

    data_r = data .* W1 .* W2
    return (data_r, W1 .* W2)
end

# function apply_rw(data::Vector{<:Array{Float64}}, W::Vector{Matrix{Float64}})
    # if length(W) != length(data)
        # error("Lenghts must match")
    # end
    # nc = size.(data, 1)
# 
    # rw1 = [W[k][1, 1:nc[k]] for k=1:length(W)]
    # rw2 = [W[k][2, 1:nc[k]] for k=1:length(W)]
    # rw = [rw1[k] .* rw2[k] for k=1:length(W)]
    # data_r = [data[k] .* rw[k] for k=1:length(data)]
    # return (data_r, rw)
# end
# rw with gaps
function apply_rw(data::Array{Float64}, W::Vector{Matrix{Float64}}, vcfg::Vector{Int64})
    idm = reshape(vcfg, :, length(W))
    nc = size.(W, 2)
    Wcat = hcat(W...)

    W1 = Wcat[1, vcfg]
    W2 = Wcat[2, vcfg]
    rw1 = [W[k][1, 1:nc[k]] for k=1:length(W)]
    rw2 = [W[k][2, 1:nc[k]] for k=1:length(W)]
    rw = [rw1[k] .* rw2[k] for k =1:length(W)]
    data_r = [data[k] .* rw[k][vcfg[k]] for k=1:length(data)]
    return (data_r, rw)
end


function plat_av(obs::Vector{uwreal}, plat::Vector{Int64}; wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)
    isnothing(wpm) ? uwerr.(obs) : [uwerr(obs[k], wpm) for k in eachindex(obs)]
    w = 1 ./ err.(obs[plat[1]:plat[2]]).^2
    av = sum(w .* obs[plat[1]:plat[2]]) / sum(w)
    return av
end

function t0_guess(t::Vector{Float64}, Ysl::Array{Float64, 3}, plat::Vector{Int64}, L::Int64)
    t2E_ax = t.^2 .* mean(mean(Ysl[:, plat[1]:plat[2], :], dims=2), dims=1)[1, 1, :] / L^3
    #look for values t2E: t2E > 0.3 and 0.0 < t2E_ax < 0.3
    if !(any(t2E_ax .> 0.3) && any((t2E_ax .< 0.3) .& (t2E_ax .> 0.0))) 
        error("Error finding solutions. Check volume normalization!")
    end
    t0_aux = minimum(abs.(t2E_ax .- 0.3))
    nt0 = findfirst(x-> abs(x - 0.3) == t0_aux, t2E_ax) #index closest value to t0
    return nt0
end

function get_model(x, p, n)
    s = 0.0
    for k = 1:n
        s = s .+ p[k] .* x.^(k-1)
    end
    return s
end