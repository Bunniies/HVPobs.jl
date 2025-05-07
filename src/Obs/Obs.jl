module Obs
    
using ..Data, ..Fits
using ADerrors, LaTeXStrings, PyPlot
using Statistics

@doc raw"""
    Corr(a::Vector{uwreal}, id::String, gamma::String)

The struct Corr stores data, ensemble id and gamma structure for
two and three-point correlation functions. 
"""
struct Corr
    obs::Vector{uwreal}
    id::String
    gamma::String

    function Corr(a::Vector{uwreal}, cd::CData)
       
        return new(a, cd.id, cd.gamma)
    end
    Corr(a::Vector{uwreal}, id::String, gamma::String) = new(a, id, gamma) 
end
function Base.show(io::IO, corr::Corr)
    println(io, "Correlator")
    println(io, " - Ensemble ID: ", corr.id)
    println(io, " - Gamma:       ", corr.gamma)
end
export Corr

@doc raw"""

"""
struct Window
    func::Function
    function Window(str::String)
        delta = 0.15

        if str == "SD" 
            d = 0.4
            @. funcsd(x0) =  1 - 0.5 * (1 + tanh((x0-d)/delta))
            return new(funcsd)
        elseif str == "ID"
            d1 = 0.4
            d2 = 1.0
            @. funcid(x0) = 0.5 * (1 + tanh((x0-d1)/delta)) - 0.5 * (1 + tanh((x0-d2)/delta))
            return new(funcid)
        elseif  str == "LD"
            d = 1.0
            @. funcld(x0) = 0.5 * (1 + tanh((x0-d)/delta))
            return new(funcld)
        elseif str == "ILD" # intermediate and long distance
            d = 0.4
            @. funcild(x0) = 0.5 * (1 + tanh((x0-d)/delta))
            return new(funcild)
        elseif str == "SID" # short and intermediate distance
            d = 2.5
            @. funcsid(x0) =  1 - 0.5 * (1 + tanh((x0-d)/delta))
            return new(funcsid)
        else
            error("Window $(str) not defined.")
        end
    end
end
function (a::Window)(x)
    return a.func(x)
end
export Window

include("ObsPrimary.jl")
export corr_obs, comp_t0, comp_fvc

include("ObsTools.jl")
export plat_av, frwd_bckwrd_symm!, frwd_bckwrd_antisymm!

include("ObsSpectrum.jl")
export meff, mpcac, dec_const

include("ObsImprovement.jl")
export improve_corr_vkvk!, improve_corr_vkvk_cons!
export ZV, cv_loc, cv_cons, bv, bv_bar, ca, Za_l_sub, ZP, cv_pert, ba_minus_bp
export ZV_set2, cv_loc_set2, cv_cons_set2, bv_set2, bv_bar_set2, ba_imp


end
