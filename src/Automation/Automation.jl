module Automation
    
    using ..Data
    using ..Obs

    struct EnsInfo
        id::String
        beta::Float64
        L::Int64
        kappa_l::Float64
        kappa_s::Float64
        kappa_crit::Float64
        plat_t0::Vector{Int64}
        dtr::Int64
        function EnsInfo(id::String)
            db = CLS_db[id]
            kappa_crit = CLS_kappa_crit[db["beta"]]
            return new(id, db["beta"], db["L"], db["kappa_l"], db["kappa_s"], kappa_crit, db["plat_t0"], db["dtr"])
        end
    end
    function Base.show(io::IO, a::EnsInfo)
        println(io, "Ensemble: ", a.id, " beta: ", a.beta, " L: ", a.L )
    end
    export EnsInfo

    include("DataAutomation.jl")
    export get_data, get_rw, get_t0, get_corr, get_fvc, get_mesons_data, get_mesons_corr

end