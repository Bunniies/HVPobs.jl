module Data

using DelimitedFiles, ADerrors, Statistics, MultiFloats

struct CData
    id::String
    rep_len::Vector{Int64}
    re_data::Array{Float64}
    im_data::Array{Float64}
    idm::Vector{Int64}
    gamma::String

    CData(id, rep_len, re_data, im_data, idm, gamma) = new(id, rep_len, re_data, im_data, idm, gamma)
end
function Base.show(io::IO, cd::CData)
    println(io, "CData")
    println(io, " - Ensemble ID:  ", cd.id)
    println(io, " - Gamma:        ", cd.gamma)
    println(io, " - Replicas:     ", length(cd.rep_len))
    println(io, " - Measurements: ", length(cd.idm))
end
export CData


struct YData
    vtr::Vector{Int32}
    t::Vector{Float64}
    obs::Array{Float64, 3}
    id::String
    YData(vtr, t, obs, id) = new(vtr, t, obs, id)
end
export YData

struct FVCData
    id::String
    n2_max::Int64
    bin::Int64
    data::Array{Float64, 3}
    FVCData(id, n2_max, bin, data) = new(id, n2_max, bin, data)
end
function Base.show(io::IO, cd::FVCData)
    println(io, "FVCata")
    println(io, " - Ensemble ID:  ", cd.id)
    println(io, " - n2_max:       ", cd.n2_max)
    println(io, " - bin:          ", cd.bin)
end
export FVCData

include("DataReader.jl")
export read_hvp_data, read_ms, read_ms1, read_FVC


include("DataConst.jl")
export GAMMA, CLS_db, hc, t0, t0sqrt_ph, CLS_kappa_crit


end
