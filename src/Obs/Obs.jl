module Obs
    
using ADerrors, LaTeXStrings, PyPlot


include("ObsTools.jl")
export plat_av

include("Spectroscopy.jl")
export meff

end
