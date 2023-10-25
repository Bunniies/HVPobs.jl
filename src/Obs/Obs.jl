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
end
function Base.show(io::IO, corr::Corr)
    println(io, "Correlator")
    println(io, " - Ensemble ID: ", corr.id)
    println(io, " - Gamma:       ", corr.gamma)
end
export Corr


include("ObsPrimary.jl")
export corr_obs, comp_t0

include("ObsTools.jl")
export plat_av

include("ObsSpectrum.jl")
export meff

include("ObsImprovement.jl")
export improve_corr_vkvk!, ZV, cv_loc

end
