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

    
    return CData(id, rep_len, re_data, im_data, idm, gamma)    
end

function read_mesons_data(path::String, id::Union{String,Nothing}=nothing)

    nn = basename(path)
    idx = occursin.(GAMMA, nn)
    gamma = GAMMA[idx][1]
    if !(gamma in GAMMA)
        error("Gamma structure $(gamma) not supported.")
    end

    f = readdlm(path, '\t',  skipstart=1)
    #return f
    header = filter(x-> typeof(x)<:AbstractString && occursin("nb", x), f)

    rep_info = last.(split.(header))
    rep_len = parse.(Int64, rep_info)
    delim = findall(x-> typeof(x)<:AbstractString && occursin("#", x), f)
    delim_line1 = delim[1:2:end]
    delim_line2 = delim[2:2:end]
    
    tvals = Int64((size(f)[1] - size(header)[1] - size(delim)[1]) / sum(rep_len))
    if tvals%2 != 0
        error("Temporal dimension T has to be an even number \n T = $(tvals)" )
    end
  

    re_data = Array{Float64}(undef, sum(rep_len), tvals)
    im_data = similar(re_data)
    
    idm = Vector{Int64}(undef, sum(rep_len))
    rep_id = Vector{Int64}(undef, sum(rep_len))

    for k in eachindex(delim_line1)

        rep_id[k], idm[k] = map(eachmatch(r"\d+", f[delim_line1[k]])) do m
            parse(Float64, m.match)
        end

        idx = delim_line2[k].I[1]
        try
            data = split.(f[idx+1:idx+tvals])
            re_data[k,:] = parse.(Float64, getindex.(data, 2))
            im_data[k,:] = parse.(Float64, getindex.(data, 3))
        catch
            re_data[k,:] = f[idx+1:idx+tvals, 2]
            im_data[k,:] = f[idx+1:idx+tvals, 3]
        end

    end
    
    if id == "E300"
        re_data = re_data[1:1137, :]
        im_data = im_data[1:1137, :]
        idm = idm[1:1137]
        rep_len = [1137]
    end
    
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
read_rw_openQCD2(path::String; print_info::Bool=false)

This function reads the reweighting factors generated with openQCD version 2.#.
The flag print_info if set to true print additional information for debugging
"""
function read_rw_openQCD2(path::String; print_info::Bool=false)
    
    data = open(path, "r")
    nrw = read(data, Int32) 
    nrw = Int32(nrw / 2)
    
    nfct = Array{Int32}(undef, nrw)
    read!(data, nfct)       
    
    nsrc = Array{Int32}(undef, nrw)
    read!(data, nsrc)        
    null = read(data, Int32)  
    if null !== Int32(0)   
        error("In openQCD 2.0 this Int32 should be a zero.")
    end

    data_array = Array{Array{Float64}}(undef, nrw)
    [data_array[k] = Array{Float64}(undef, 0) for k in 1:nrw]
    vcfg = Vector{Int32}(undef, 0) 
    while !eof(data)

        push!(vcfg, read(data, Int32)) 
        if print_info
            println("\n ######## cnfg: ", vcfg[end])
        end

        for k in 1:nrw
            read_array_rwf_dat_openQCD2(data)
            tmp_rw, n = read_array_rwf_dat_openQCD2(data)

            tmp_nfct=1.0
            for j in 1:n[1]
                tmp_nfct *= mean((exp.(.-tmp_rw[j])))
            end
            push!(data_array[k], tmp_nfct)
        end
    end

    # return data_array
    return permutedims(hcat(data_array...), (2,1))
end


function read_array_rwf_dat_openQCD2(data::IOStream; print_info::Bool=false)
    
    d = read(data, Int32) 
    n = Array{Int32}(undef, d)
    read!(data, n)   
    size = read(data, Int32) 
    m = prod(n)
    
    if print_info
        println("d: ", d)
        println("n: ", n)
        println("size: ", size)
        println("m: ", m)
    end

    if size == 4
        types = Int32
    elseif size == 8
        types = Float64
    elseif size == 16
        types = Float64x2
    else
        error("No type with size=$(size) supported")
    end

    tmp_data = Array{types}(undef, m)
    read!(data, tmp_data)

    res = parse_array_openQCD2(d, n, tmp_data, quadprec=true) 

    return res, n
end

function parse_array_openQCD2(d, n, dat; quadprec=true)
    
    if d != 2
        error("dat must be a two-dimensional array")
    end
    res = Vector{Vector{Float64}}(undef, 0)
    
    for k in range(1,n[1])
        tmp = dat[(k-1)*n[2]+1:k*n[2]]
        if quadprec
            tmp2 = Vector{Float64}(undef, 0)
            for j in range(start=1,step=2,stop=length(tmp))
                push!(tmp2, tmp[j])
            end
            push!(res, tmp2)
        else
            push!(res, tmp)
        end
    end

    return res
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

    if v == "2.0"
        return read_rw_openQCD2(path)
    end
    r_data = read_rw(path, v=v)
    nrw = length(r_data)
    ncnfg = size(r_data[1])[3]
    W = zeros(Float64, nrw, ncnfg)
    [W[k, :] = prod(mean(exp.(.-r_data[k]), dims=2), dims=1) for k = 1:nrw]
    return W
end

function read_FVC(path::String)
    
    id = match(r"[A-Z]+[0-9]{3}", basename(path)).match
    f = readdlm(path, '\t', '\n' )
    header = filter(x->  typeof(x)<:AbstractString &&occursin("#", x), f )

    n2_max = parse(Int64, last(split(filter(x->occursin("n2_max",x), header)[1])))
    info = split(filter(x->occursin("n_t_lat",x), header)[1])
    Thalf = Int64(parse(Int64, info[end-2])/2)
    bin = parse(Int64, info[end-3])
   
    # data = setdiff(f, header)
    data = filter(x -> typeof(x)==Float64, f)

    jkdata = Array{Float64}(undef, n2_max, Thalf, bin )

    for n2 in 1:n2_max
        for t in 1:Thalf
            for b in 1:bin
               jkdata[n2, t, b] = data[b + bin*(t-1) + bin*Thalf*(n2-1)]
            end
        end
    end

    return FVCData(id, n2_max, bin, jkdata)

end
