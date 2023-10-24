module HVPobs

include("Data/Data.jl")

using .Data 
export CData, Corr, YData
export GAMMA
export read_hvp_data, read_ms, read_ms1

include("Fits/Fits.jl")

using .Fits
export FitRes
export get_pvalue
export fit_routine

include("Obs/Obs.jl")

using .Obs
export meff 
export corr_obs, comp_t0


end # module
