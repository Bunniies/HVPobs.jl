
function get_data(path::String, ens::String, fl::String, g::String)

    if !(fl in ["light", "pion", "strange", "charm", "charm_plus"])
        error("Flavour fl not found. \n Choose fl from: light, pion, strange, charm, charm_plus ")
    end
    if !(g in GAMMA)
        error("Gamma structure \"$(g)\" not found in $(GAMMA)")
    end
    p = joinpath(path, ens, fl, "raw_data")

    data = filter(x-> last(split(basename(x), "_")) == g, readdir(p, join=true))
    
    if isempty(data)
        error("No data found for ensemble $(ens) with flavour \"$(fl)\" and gamma structure  \"$(g)\" ")
    end

    return read_hvp_data(data[1], ens)
end

function get_rw(path::String, ens::String)
    p = joinpath(path, ens)
    rep = filter(x->occursin("ms1.dat", x), readdir(p, join=true))
    if ens == "J500"
        return [read_ms1(rep[1]), read_ms1(rep[2], v="1.4")]
    end
    if ens == "J501"
        return [read_ms1(rep[1]), read_ms1(rep[2], v="1.4") ]
    end
    if ens == "E250"
        return read_ms1(rep[1], v="1.6")
    end
    if length(rep)!=0
        try
            length(rep) == 1 ? (return read_ms1(rep[1])) : (return read_ms1.(rep)) 
        catch
            length(rep) == 1 && length(rep)!=0 ? (return read_ms1(rep[1], v="1.4")) : (return read_ms1.(rep, v="1.4")) 
        end
    else
        error("ms1.dat file not found for ensemble ", ens, " in path ", p)
    end
end

function get_t0(path::String, ens::String, dtr::Int64)
    p = joinpath(path, ens)
    rep = filter(x->occursin("ms.dat", x), readdir(p, join=true))
    if length(rep)!=0
        length(rep) == 1 ? (return read_ms(rep[1], dtr=dtr, id=ens)) : (return read_ms.(rep, dtr=dtr, id=ens)) 
    else
        error("ms.dat file not found for ensemble ", ens, " in path ", p)
    end
end

function get_corr(path::String, ens::String, fl::String, g::String, path_rw::Union{String, Nothing}=nothing)

    cdata = get_data(path, ens, fl, g)
    rw = isnothing(path_rw) ? nothing : get_rw(path_rw, ens)
    return corr_obs(cdata, real=true, rw=rw)
end