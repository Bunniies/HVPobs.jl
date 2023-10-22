function corr_obs(cd::CData; real::Bool=true, rw::Union{Array{Float64,2}, Nothing}=nothing, L::Int64=1, nms::Union{Int64, Nothing}=nothing)
    
    real ? data = cd.re_data ./ L^3 : data = cd.im_data ./ L^3

    cnfg = size(data, 1)
    tvals = size(data, 2)
    nms = isnothing(nms) ?  sum(cd.rep_len) : nms


    if isnothing(rw)
        obs = [uwreal(data[:,t], cd.id, cd.rep_len, cd.idm, nms) for t in 1:tvals]
    else
        # code reweighting 
    end

    return Corr(obs, cd)
end