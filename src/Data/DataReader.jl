function read_hvp_data(path::String, id::Union{String,Nothing}=nothing)

    nn = basename(path)
    gamma = split(nn, "_")[end]
    if !(gamma in GAMMA)
        error("Gamma structure $(gamma) not supported.")
    end

    f = readdlm(path, '\t', '\n', skipstart=1)
    
    header = filter(x-> typeof(x)<:AbstractString && occursin("nb", x), f)

    rep_info = last.(split.(header))
    rep_len = parse.(Int64, rep_info)
    delim = findall(x-> typeof(x)<:AbstractString && occursin("#", x), f)

    tvals = Int64((size(f)[1] - size(header)[1] - size(delim)[1]) / sum(rep_len))
    if tvals%2 != 0
        error("Temporal dimension T has to be an even number \n T = $(tvals)" )
    end

    re_data = Array{Float64}(undef, sum(rep_len), tvals)
    im_data = similar(re_data)
    
    idm = Vector{Int64}(undef, sum(rep_len))
    rep_id = Vector{Int64}(undef, sum(rep_len))

    for k in eachindex(delim)

        rep_id[k], idm[k] = map(eachmatch(r"\d+", f[delim[k]])) do m
            parse(Float64, m.match)
        end

        idx = delim[k].I[1]
        re_data[k,:] = f[idx+1:idx+tvals, 2]
        im_data[k,:] = f[idx+1:idx+tvals, 3]
    end

    if length(rep_len) != 1
        idm_sum = [fill( (k-1)*sum(rep_len[1:k-1]), rep_len[k]) for k in eachindex(rep_len)]
    end
    # for k in eachindex(rep_len)
        # fill( (k-1)*sum(rep_len[1:k-1]), rep_len[k])
    # end
    idm .+= vcat(idm_sum...)
    return CData(id, rep_len, re_data, im_data, idm, gamma)    
end

