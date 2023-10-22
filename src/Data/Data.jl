module Data

using DelimitedFiles, ADerrors

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

struct Corr
    obs::Vector{uwreal}
    id::String
    gamma::String

    function Corr(a::Vector{uwreal}, cd::CData)
       
        return new(a, cd.id, cd.gamma)
    end
end
function Base.show(io::IO, corr::Corr)
    println(io, "Correlator")
    println(io, " - Ensemble ID: ", corr.id)
    println(io, " - Gamma:       ", corr.gamma)
end
export Corr


include("DataReader.jl")
export read_hvp_data

include("DataObs.jl")
export corr_obs

include("DataConst.jl")
export GAMMA


end
