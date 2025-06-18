module Data

using DelimitedFiles, ADerrors, Statistics, MultiFloats, PyCall, OrderedCollections

mutable struct CData
    id::String
    rep_len::OrderedDict{String,Int64}
    re_data::Array{Float64}
    im_data::Array{Float64}
    idm::Vector{Int64}
    gamma::String
    nms::Int64
    replicatot::OrderedDict{String,Int64}

    kappas::Union{Nothing, Vector{Float64}}
    theta1::Union{Nothing, Vector{Float64}}
    theta2::Union{Nothing, Vector{Float64}}

    CData(id, rep_len, re_data, im_data, idm, gamma) = new(id, rep_len, re_data, im_data, idm, gamma, CLS_CNFG[id]["nms"], CLS_CNFG[id]["repLen"] )
    CData(id, rep_len, re_data, im_data, idm, gamma, kappas, theta1, theta2) = new(id, rep_len, re_data, im_data, idm, gamma, CLS_CNFG[id]["nms"], CLS_CNFG[id]["repLen"], kappas, theta1, theta2 )

end
function Base.show(io::IO, cd::CData)
    println(io, "CData")
    println(io, " - Ensemble ID:  ", cd.id)
    println(io, " - Gamma:        ", cd.gamma)
    println(io, " - Replicas:     ", length(cd.rep_len))
    println(io, " - Measurements: ", length(cd.idm))
end
export CData


mutable struct YData
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
export read_hvp_data, read_mesons_data, read_ms, read_ms1, read_FVC, read_tree_level_v33, read_tree_level_v3sig03
export read_disconnected_from_Marcos_data, read_disconnected_from_npz, read_rwf_strange, read_kappa_charm_all_config, get_kappa_values
export read_Bphysics_data
export concat_data!, truncate_data!
export  print_uwreal

include("DataConst.jl")
export GAMMA, CLS_db, CLS_kappa_crit, Zvc_l


end
