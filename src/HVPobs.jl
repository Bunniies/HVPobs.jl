module HVPobs

include("LinAlg/LinAlg.jl")

using .LinAlg

include("Data/Data.jl")

using .Data 
export CData, Corr, YData
export GAMMA, CLS_db, hc, t0, t0sqrt_ph, CLS_kappa_crit
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
export improve_corr_vkvk!, ZV, cv_loc, cv_cons, bv, bv_bar, improve_corr_vkvk_cons!
export frwd_bckwrd_symm!

include("Automation/Automation.jl")

using .Automation
export EnsInfo
export get_data, get_rw, get_t0, get_corr


end # module
