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

    # if length(rep_len) != 1
        # idm_sum = [fill((k-1)*sum(rep_len[1:k-1]), rep_len[k]) for k in eachindex(rep_len)]
        # idm .+= vcat(idm_sum...)
    # end
    
    return CData(id, rep_len, re_data, im_data, idm, gamma)    
end

@doc raw"""
    read_ms(path::String; id::Union{String, Nothing}=nothing, dtr::Int64=1, obs::String="Y")  

Reads openQCD ms dat files at a given path. This method return YData: 
    
- `t(t)`: flow time values

- `obs(icfg, x0, t)`: the time-slice sums of the densities of the observable (Wsl, Ysl or Qsl)

- `vtr`: vector that contains trajectory number

- `id`: ensemble id

`dtr` = `dtr_cnfg` / `dtr_ms`, where `dtr_cnfg` is the number of trajectories computed before saving the configuration. `dtr_ms`
is the same but applied to the ms.dat file.

Examples:
```@example    
Y = read_ms(path)
```
"""
function read_ms(path::String; id::Union{String, Nothing}=nothing, dtr::Int64=1 , obs::String="Y")  
    if isnothing(id)
        bname = basename(path)
        m = findfirst(r"[A-Z][0-9]{3}r[0-9]{3}", bname)
        id = bname[m[1:4]]
        #id = parse(Int64, bname[m[2:4]])
    end  
    data = open(path, "r")

    dn = read(data, Int32)
    nn = read(data, Int32)
    tvals = read(data, Int32)
    eps = read(data, Float64)

    fsize = filesize(path)
    datsize=4 + 3*8*(nn + 1) * tvals # measurement size of each cnfg

    ntr = Int32((fsize - 3*4 - 8) / datsize)

    if mod(ntr, dtr) != 0
        println("ntr = ", ntr)
        println("dtr = ", dtr)
        error("ntr / dtr must be exact")
    end

    vntr = Vector{Int32}(undef, div(ntr, dtr))
    
    # x0, t, cfg 
    Wsl = Array{Float64}(undef, div(ntr, dtr), tvals, nn + 1) 
    Ysl = Array{Float64}(undef, div(ntr, dtr), tvals, nn + 1)
    Qsl = Array{Float64}(undef, div(ntr, dtr), tvals, nn + 1)

    k = 0
    for itr = 1:ntr
        tmp = read(data, Int32)
        if mod(itr, dtr) == 0
            k += 1
            vntr[k] = tmp
        end

        for iobs = 1:3
            for inn = 0:nn
                tmp2 = Vector{Float64}(undef, tvals)
                read!(data, tmp2)
                if mod(itr, dtr) == 0
                    if iobs == 1
                        Wsl[k, :, inn + 1] = tmp2
                    elseif iobs == 2
                        Ysl[k, :, inn + 1] = tmp2
                    elseif iobs == 3
                        Qsl[k, :, inn + 1] = tmp2
                    end
                end
            end
        end
    end
    close(data)
    t = Float64.(0:nn) .* dn .* eps
    
    if obs == "W"
        return YData(vntr, t, Wsl, id)
    elseif obs == "Y"
        return YData(vntr, t, Ysl, id)
    elseif obs == "Q"
        return YData(vntr, t, Qsl, id)
    else
        error("obs = ", obs," is not valid. \n Please choose between W,Y,Q")
        return nothing
    end
end


function read_rw(path::String; v::String="1.2")
    data = open(path, "r")
    nrw = read(data, Int32)

    if v == "1.4" || v =="1.6"
        nfct = Array{Int32}(undef, nrw)
        read!(data, nfct)
        nfct_inheader = 1
    elseif v=="1.2"
        nfct = ones(Int32, nrw)
        nfct_inheader = 0
    else
        error("Version not supported")  
    end
    nsrc = Array{Int32}(undef, nrw)
    read!(data, nsrc)
    glob_head_size = 4 + 4*nrw*(nfct_inheader + 1)
    datsize = 4 + 2*8*sum(nsrc .* nfct)

    fsize = filesize(path)

    ncnfg = Int32((fsize - glob_head_size)/datsize)
    r_data = Array{Array{Float64}}(undef, nrw)
    
    [r_data[k] = zeros(Float64, nfct[k], nsrc[k], ncnfg) for k = 1:nrw]
    vcfg = Array{Int32}(undef, ncnfg)
    for icfg = 1:ncnfg
        vcfg[icfg] = read(data, Int32)
        for irw = 1:nrw
            for ifct = 1:nfct[irw]
                tmp = zeros(Float64, nsrc[irw])
                seek(data, position(data) + 8 * nsrc[irw])
                read!(data, tmp)
                r_data[irw][ifct, :, icfg] = tmp[:]
            end
        end
    end
    close(data)
    return r_data
end

@doc raw"""
    read_ms1(path::String; v::String="1.2")

Reads openQCD ms1 dat files at a given path. This method returns a matrix `W[irw, icfg]` that contains the reweighting factors, where
`irw` is the `rwf` index and icfg the configuration number.
The function is compatible with the output files of openQCD v=1.2, 1.4 and 1.6. Version can be specified as argument.

Examples:
```@example    
read_ms1(path)
read_ms1(path, v="1.4")
read_ms1(path, v="1.6")
```
"""
function read_ms1(path::String; v::String="1.2")
    r_data = read_rw(path, v=v)
    nrw = length(r_data)
    ncnfg = size(r_data[1])[3]
    W = zeros(Float64, nrw, ncnfg)
    [W[k, :] = prod(mean(exp.(.-r_data[k]), dims=2), dims=1) for k = 1:nrw]
    return W
end