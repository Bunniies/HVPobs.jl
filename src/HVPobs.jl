module HVPobs



include("LinAlg/LinAlg.jl")
export Rationaluw, uwerr
export abs

using .LinAlg

include("Data/Data.jl")

using .Data 
export CData, Corr, YData
export GAMMA, CLS_db, CLS_kappa_crit # hc, t0, a, , mpi_ph, Zvc_l
export read_hvp_data, read_mesons_data, read_ms, read_ms1, read_FVC, read_tree_level_v33, read_tree_level_v3sig03, read_disconnected_from_Marcos_data, read_disconnected_from_npz

include("Fits/Fits.jl")

using .Fits
export FitRes
export get_pvalue
export fit_routine
export bayesian_av

include("Obs/Obs.jl")

using .Obs
export Corr, Window
export meff, mpcac 
export corr_obs, comp_t0, comp_fvc
export improve_corr_vkvk!, ZV, cv_loc, cv_cons, bv, bv_bar, improve_corr_vkvk_cons!, cv_loc_set2, cv_cons_set2, bv_set2, ZV_set2, bv_bar_set2, ca, Za_l_sub, ZP
export plat_av, frwd_bckwrd_symm!

include("Automation/Automation.jl")

using .Automation
export EnsInfo
export get_data, get_rw, get_t0, get_corr, get_fvc, get_mesons_data, get_mesons_corr, get_data_disc

include("LMA/LMA.jl")

using .LMA
export read_eigen_eigen

end # module
