module Data

using DelimitedFiles, ADerrors, Statistics

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




include("DataReader.jl")
export read_hvp_data, read_ms, read_ms1


include("DataConst.jl")
export GAMMA


end
