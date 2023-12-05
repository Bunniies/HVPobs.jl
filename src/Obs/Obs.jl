module Obs
    
using ..Data, ..Fits
using ADerrors, LaTeXStrings, PyPlot
using Statistics

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


include("ObsPrimary.jl")
export corr_obs, comp_t0, comp_fvc

include("ObsTools.jl")
export plat_av, frwd_bckwrd_symm!

include("ObsSpectrum.jl")
export meff

include("ObsImprovement.jl")
export improve_corr_vkvk!, ZV, cv_loc, cv_cons, bv, bv_bar, improve_corr_vkvk_cons!

end
