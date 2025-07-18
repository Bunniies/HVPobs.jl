
function get_data(path::String, ens::String, fl::String, g::String)

    if !(fl in ["light", "light_LMA", "pion", "strange", "charm", "charm_plus"])
        error("Flavour $(fl) not found. \n Choose fl from: light, pion, strange, charm, charm_plus ")
    end
    if !(g in GAMMA)
        error("Gamma structure \"$(g)\" not found in $(GAMMA)")
    end

    
    if fl == "light"
        try
            p = joinpath(path, ens, "light_LMA", "raw_data")
            data = filter(x-> last(split(basename(x), "_")) == g, readdir(p, join=true))
            return read_hvp_data(data[1], ens)
        catch
            p = joinpath(path, ens, fl, "raw_data")
            data = filter(x-> last(split(basename(x), "_")) == g, readdir(p, join=true))
            return read_hvp_data(data[1], ens)
        end
    elseif  fl == "pion" 
        p = joinpath(path, ens, fl)        
    else 
        p = joinpath(path, ens, fl, "raw_data")
    end

    data = filter(x-> last(split(basename(x), "_")) == g, readdir(p, join=true))
    if isempty(data)
        error("No data found in $(p) for ensemble $(ens) with flavour \"$(fl)\" and gamma structure  \"$(g)\" ")
    end

    return read_hvp_data(data[1], ens)
end

function get_data_Bphysics(path::String, ens::String, fl::String, g::String)
    if !(g in GAMMA)
        error("Gamma structure \"$(g)\" not found in $(GAMMA)")
    end

    path = joinpath(path, fl)
    fdata = filter(x-> occursin(g, x), readdir(path, join=true))
    if isempty(fdata)
        error("No data file found for Ens $(ens), sector $(fl), gamma $(g)")
    end
    return read_Bphysics_data(fdata[1], ens)
end


function get_data_disc(path::String, ens::String, fl::String)
    if fl ∉ ["08", "0c", "80", "88", "8c", "c0", "c8", "cc"]
        error("Unrecognised flavour structure $(fl): choose from [08, 0c, 80, 88, 8c, c0, c8, cc]")
    end

    p = joinpath(path, "disc", ens, "disc/2pt")

    data = filter(x-> occursin("$(fl).npz", basename(x)), readdir(p, join=true))

    if isempty(data)
        error("No data found in $(p) for ensemble $(ens) with flavour \"$(fl)\" ")
    end
    return read_disconnected_from_npz(data[1], ens)    
end

function get_mesons_data(path::String, ens::String, fl::String, g::String )
    if !(fl in ["ll", "ls", "ss"])
        error("Flavour $(fl) not found. Choose fl from: ll, ls, ss")
    end
    if !(g in GAMMA)
        error("Gamma structure \"$(g)\" not found in $(GAMMA)")
    end

    p = joinpath(path, ens)

    data = filter(x-> occursin(fl, basename(x)) && occursin(g, basename(x)), readdir(p,join=true))
    if isempty(data)
        error("No data found for ensemble $(ens) with flavour \"$(fl)\" and gamma structure  \"$(g)\" ")
    end

    return read_mesons_data(data[1], ens)
end

function get_rw(path::String, ens::String; v::String="1.2")
    rep = filter(x->occursin(ens, x), readdir(path, join=true))
    if ens == "J500"
        return [read_ms1(rep[1]), read_ms1(rep[2], v="1.4"), read_ms1(rep[3], v="2.0") ]
    end
    if ens == "J501"
        return [read_ms1(rep[1]), read_ms1(rep[2], v="1.4"), read_ms1(rep[3], v="2.0") ]
    end
    if ens == "D450"
        return [read_ms1(rep[1], v="2.0"), read_ms1(rep[2], v="2.0")]
    end
    if ens == "E300"
        return [read_ms1(rep[1], v="2.0"), read_ms1(rep[2], v="2.0"), read_ms1(rep[3], v="2.0")] 
    end
    # if ens == "N202" # TO USE ONLY IN B PHYSICS PROJECT!
    #     return [read_ms1(rep[1]), read_ms1(rep[2], v="1.4")]
    # end
    if length(rep)!=0
        try
            length(rep) == 1 ? (return read_ms1(rep[1], v=v)) : (return read_ms1.(rep, v=v)) 
        catch
            try
                length(rep) == 1 && length(rep)!=0 ? (return read_ms1(rep[1], v="1.4")) : (return read_ms1.(rep, v="1.4")) 
            catch
                length(rep) == 1 && length(rep)!=0 ? (return read_ms1(rep[1], v="2.0")) : (return read_ms1.(rep, v="2.0")) 
            end
        end
    else
        error("ms1.dat file not found for ensemble ", ens, " in path ", p)
    end
end

function get_corr(path::String, ens::EnsInfo, fl::String, g::String; path_rw::Union{String, Nothing}=nothing, L::Int64=1, frw_bcwd::Bool=false)

    cdata = get_data(path, ens.id, fl, g)
    rw = isnothing(path_rw) ? nothing : get_rw(path_rw, ens.id)
    corr = corr_obs(cdata, real=true, rw=rw, L=L)
    if frw_bcwd
        frwd_bckwrd_symm!(corr)
    end

    return corr 
end

function get_Bphysics_corr(path::String, ens::EnsInfo, fl::String, g::String; path_rw::Union{String, Nothing}=nothing, L::Int64=1)
    
    cdata = get_data_Bphysics(path, ens.id, fl, g)
    rw = isnothing(path_rw) ? nothing : get_rw(path_rw, ens.id)
    corr = corr_obs(cdata, real=true, rw=rw, L=L)

    return corr
end

function get_corr_disc(path::String, ens::EnsInfo, fl::String ; path_rw::Union{String, Nothing}=nothing, L::Int64=1, frw_bcwd::Bool=false)
    if fl ∉ ["08", "0c", "80", "88", "8c", "c0", "c8", "cc"]
        error("Unrecognised flavour structure $(fl): choose from [08, 0c, 80, 88, 8c, c0, c8, cc]")
    end

    cdata = get_data_disc(path, ens.id, fl)
    KEYS = collect(keys(cdata))
    rw = isnothing(path_rw) ? nothing : get_rw(path_rw, ens.id)
    
    corr_dict  = OrderedDict{String,Corr}()
    for kk in KEYS
        corr = corr_obs(cdata[kk], real=true, rw=rw, L=L)
        if frw_bcwd
            frwd_bckwrd_symm!(corr)
        end
        corr_dict[kk] = corr
    end
    return corr_dict
end

function get_mesons_corr(path::String, ens::EnsInfo, fl::String, g::String; path_rw::Union{String, Nothing}=nothing, L::Int64=1, frw_bcwd::Bool=false)

    cdata = get_mesons_data(path, ens.id, fl, g)
    if ens.id == "A653"
        println("Ensemble A653: rwf are taken from ens A654 as the correct rwf are not currently available")
        rw = isnothing(path_rw) ? nothing : get_rw(path_rw, "A654")
    else
        rw = isnothing(path_rw) ? nothing : get_rw(path_rw, ens.id)
    end
    corr = corr_obs(cdata, real=true, rw=rw, L=L)
    if frw_bcwd
        frwd_bckwrd_symm!(corr)
    end

    return corr 
end

function get_t0(path::String, ens::EnsInfo; pl::Bool=false, path_rw::Union{String, Nothing}=nothing, wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)

    data = read_t0(path, ens.id, ens.dtr)
    if ens.id == "A653"
        println("Ensemble A653: rwf are taken from ens A654 as the correct rwf are not currently available")
        rw = isnothing(path_rw) ? nothing : get_rw(path_rw, "A654")
    else
        rw = isnothing(path_rw) ? nothing : get_rw(path_rw, ens.id)
    end
    if ens.id == "F300"
        truncate_data!(data[1], 105)
    end
    t0_res = comp_t0(data, ens.plat_t0, L=ens.L, pl=pl, rw=rw, info=false, wpm=wpm)
    return t0_res
end


function read_t0(path::String, ens::String, dtr::Int64)
    # p = joinpath(path, ens)
    rep = filter(x->occursin(ens, x), readdir(path, join=true))
    if length(rep)!=0
        length(rep) == 1 ? (return read_ms(rep[1], dtr=dtr, id=ens)) : (return read_ms.(rep, dtr=dtr, id=ens)) 
    else
        error("ms.dat file not found for ensemble ", ens, " in path ", path)
    end
end

function get_fvc(path::String, ens::String)
    rep = filter(x-> occursin("corr_blat_gsl", x) && occursin(ens, x), readdir(path, join=true))[1]
    return comp_fvc(rep)
end


