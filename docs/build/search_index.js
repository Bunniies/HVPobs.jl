var documenterSearchIndex = {"docs":
[{"location":"fits/index_fits.html#Fits-Module","page":"Fits","title":"Fits Module","text":"","category":"section"},{"location":"fits/index_fits.html","page":"Fits","title":"Fits","text":"FitRes\nfit_routine\nget_pvalue\nbayesian_av","category":"page"},{"location":"fits/index_fits.html#HVPobs.Fits.FitRes","page":"Fits","title":"HVPobs.Fits.FitRes","text":"@doc raw     FitRes structure.\n\nThis structure is returned by the fit_routine function and it contains useful fit informations such as:       - dof       - param       -chi2       -chi2exp       -pval  \n\n\n\n\n\n","category":"type"},{"location":"fits/index_fits.html#HVPobs.Fits.fit_routine","page":"Fits","title":"HVPobs.Fits.fit_routine","text":"fit_routine(model::Function, xdata::Array{<:Real}, ydata::Array{uwreal}, param::Int64=3; wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)\nfit_routine(model::Vector{Function}, xdata::Vector{Array{Float64, N}} where N, ydata::Vector{Array{uwreal, N}} where N, param::Int64; wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing, correlated_fit::Bool=false, pval::Bool=false)\n\nGiven a model function with a number param of parameters,  an array of xdata and an array of ydata of type uwreal, this function performs the fit with the given model and returns a FitRes object.\n\nIf an array of model is passed instead, together with a vector of array of xdata and ydata, the function performs a global fit by applying each  model element-wise to each set of data.\n\nThis function automatically performs error propagation from ydata to the fit parameters returned in FitRes.\n\nThe flag pval allows to compute the fit pvalues, if set to true.\n\nThe flag wpm is an optional array of Float64 of lenght 4 used by ADerrors package for error analysis. The first three paramenters specify the criteria to determine the summation windows:\n\nvp[1]: The autocorrelation function is summed up to t = round(vp1).\nvp[2]: The sumation window is determined using U. Wolff poposal with S_tau = wpm2\nvp[3]: The autocorrelation function Gamma(t) is summed up a point where its error deltaGamma(t) is a factor vp[3] times larger than the signal.\n\nAn additional fourth parameter vp[4], tells ADerrors to add a tail to the error with tau_exp = wpm4. Negative values of wpm[1:4] are ignored and only one component of wpm[1:3] needs to be positive. If the flag covaris set to true, fit_routine takes into account covariances between x and y for each data point.\n\n@. model(x,p) = p[1] + p[2] * exp(-(p[3]-p[1])*x)\n@. model2(x,p) = p[1] + p[2] * x[:, 1] + (p[3] + p[4] * x[:, 1]) * x[:, 2]\n\n# single fit\nfit_routine(model, xdata, ydata, param=3)\n# global fit\nfit_routine([model, model2], [xdata, xdata2], [ydata, ydata2], param=3)\n\n\n\n\n\n","category":"function"},{"location":"fits/index_fits.html#HVPobs.Fits.get_pvalue","page":"Fits","title":"HVPobs.Fits.get_pvalue","text":"pvalue(chisq::Function,\nchi2::Float64,\nxp::Vector{Float64}, \ndata::Vector{uwreal},\nwpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing} = Dict{Int64,Vector{Float64}}();\nW::Union{Vector{Float64},Array{Float64,2}} = Vector{Float64}(),\nnmc::Int64 = 5000)\n\nComputes the p-value of a previously done fit, using as input the \\chi^2 observed from the fit, the fit parameters and the fitted data.  The p-value for a given \\chi^2 is the probability of, given the data you have, finding such a \\chi^2 or worse from a fit, and still have the data well described by the fit function. nmc is the number of MC samples used to estimate the p-value integral, default is 5000. By now it only works with a vector for weights (containing the diagonal of W)\n\nfunction fit_defs(f::Function,x,W)\n    chisq(p,d)=sum((d-f(x,p)).^2 .*W)\n    return chisq\nend\n\n@.fun(x,p) = p[1] + p[2] * x\nchisq = fit_defs(fun, x, 1.0 ./ err.(y) .^ 2)\nfit = curve_fit(fun, x, value.(y), 1.0 ./ err.(y) .^ 2, [0.5, 0.5])\n(up, chiexp) = fit_error(chisq, coef(fit), y)\n\nwpm = Dict{Int64,Vector{Float64}}()\nwpm[1] = [-1.0,-1.0,4-0,-1.0]\nQ = pvalue(chisq, chi2, value.(up), y, wpm; W = 1.0 ./ err.(y) .^ 2, nmc=10000)\n#Q = pvalue(chisq, chi2, value.(up), y; W = 1.0 ./ err.(y) .^ 2, nmc=10000)\n#Q = pvalue(chisq, chi2, value.(up), y)\n\n\n\n\n\n","category":"function"},{"location":"fits/index_fits.html#HVPobs.Fits.bayesian_av","page":"Fits","title":"HVPobs.Fits.bayesian_av","text":"bayesian_av(fun::Function, y::Array{uwreal}, tmin_array::Array{Int64}, tmax_array::Array{Int64}, k::Int64, pl::Bool, data::Bool; wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)\n\nbayesian_av(fun1::Function, fun2::Function, y::Array{uwreal}, tmin_array::Array{Int64}, tmax_array::Array{Int64}, k1::Int64, k2::Int64, pl::Bool, data::Bool; wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)\n\nbayesian_av(fun::Array{Function}, y::Array{uwreal}, tmin_array::Array{Int64}, tmax_array::Array{Int64}, k::Array{Int64}, pl::Bool, data::Bool; wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)\n\nComputes bayesian average of data. For a given fit function, it explores choices of fit intervals, assigning each of them a weight. The function saves the first fit parameter of your function, and then it does the weighted average of it and assigns a systematic. See https://arxiv.org/abs/2008.01069 \n\nThe function takes as input the fit intervals to explore. \n\ntmin_array is an array of integers with the lower bounds on the fit intervals to explore, ***ordered from lower to higher***.\n\ntmax_array is an array of integers with the upper bounds on the fit intervals to explore, ***ordered from lower to higher***.\n\nk is the number of parameters of the fit function to use.\n\nYou can also use as input two fit functions, and two values of k, one for each function. Then, for each fit interval choice, the function explores the two fit functions. This means that for each fit interval choice you get two results: one for the first fit funcction, and another for the second. You can also use a vector of functions and a vector of k (numer of parameters of each funtion) to apply the bayesian averaging method to multiple functions.\n\nThe method returns two objects: first, the weighted average as an uwreal object, with mean value and statistichal error. The second object returned is the systematic error coming from the fit interval variation. If data is true, then returns 4 objects: weighted average, systematic error, a vector with the results of the fit for each fit interval choice, and a vector with the weights associated to each fit.\n\n@.fun(x,p) = p[1] * x ^0\nk = 1\ntmin_array = [10,11,12,13,14,15]\ntmax_array = [80,81,82,83,84,85]\n(average, systematics, data, weights) = bayesian_av(fun,x,tmin_array,tmax_array,k,pl=true,data=true)\n\n@.fun1(x,p) = p[1] * x ^0\n@.fun2(x,p) = p[1] + p[2] * exp( - p[3] * (x))\nk1 = 1\nk2 = 3\ntmin_array = [10,11,12,13,14,15]\ntmax_array = [80,81,82,83,84,85]\n(average, systematics) = bayesian_av(fun1,fun2,x,tmin_array,tmax_array,k1,k2)\n\n@.fun1(x,p) = p[1] * x ^0\n@.fun2(x,p) = p[1] + p[2] * exp( - p[3] * (x))\nk1 = 1\nk2 = 3\ntmin_array = [10,11,12,13,14,15]\ntmax_array = [80,81,82,83,84,85]\n(average, systematics) = bayesian_av([fun1,fun2],x,tmin_array,tmax_array,[k1,k2])\n\n\n\n\n\n","category":"function"},{"location":"data/index_data.html#Data-Module","page":"Data","title":"Data Module","text":"","category":"section"},{"location":"data/index_data.html","page":"Data","title":"Data","text":"read_hvp_data\nread_mesons_data\nread_ms\nread_ms1","category":"page"},{"location":"data/index_data.html#HVPobs.Data.read_hvp_data","page":"Data","title":"HVPobs.Data.read_hvp_data","text":"@doc raw     readhvpdata(path::String, id::String)\n\nThis function reads the HVP data used in the g-2 analysis. It takes as input a path to the data and en ensemble id  and it returns a CData structure.\n\n@example cdata = read_hvp_data(path, \"H101\")`\n\n\n\n\n\n","category":"function"},{"location":"data/index_data.html#HVPobs.Data.read_mesons_data","page":"Data","title":"HVPobs.Data.read_mesons_data","text":"@doc raw     readmesondata(path::String, id::String)\n\nThis function reads the 2-point function data in the format used for spectrum computations. It takes as input a path to the data and en ensemble id  and it returns a CData structure.\n\n@example cdata = read_mesons_data(path, \"H101\")`\n\n\n\n\n\n","category":"function"},{"location":"data/index_data.html#HVPobs.Data.read_ms","page":"Data","title":"HVPobs.Data.read_ms","text":"read_ms(path::String; id::Union{String, Nothing}=nothing, dtr::Int64=1, obs::String=\"Y\")\n\nReads openQCD ms dat files at a given path. This method return YData: \n\nt(t): flow time values\nobs(icfg, x0, t): the time-slice sums of the densities of the observable (Wsl, Ysl or Qsl)\nvtr: vector that contains trajectory number\nid: ensemble id\n\ndtr = dtr_cnfg / dtr_ms, where dtr_cnfg is the number of trajectories computed before saving the configuration. dtr_ms is the same but applied to the ms.dat file.\n\nExamples:\n\nY = read_ms(path)\n\n\n\n\n\n","category":"function"},{"location":"data/index_data.html#HVPobs.Data.read_ms1","page":"Data","title":"HVPobs.Data.read_ms1","text":"read_ms1(path::String; v::String=\"1.2\")\n\nReads openQCD ms1 dat files at a given path. This method returns a matrix W[irw, icfg] that contains the reweighting factors, where irw is the rwf index and icfg the configuration number. The function is compatible with the output files of openQCD v=1.2, 1.4 and 1.6. Version can be specified as argument.\n\nExamples:\n\nread_ms1(path)\nread_ms1(path, v=\"1.4\")\nread_ms1(path, v=\"1.6\")\nread_ms1(path, v=\"2.0\")\n\n\n\n\n\n\n","category":"function"},{"location":"obs/index_obs.html#Obs-Module","page":"Obs","title":"Obs Module","text":"","category":"section"},{"location":"obs/index_obs.html#ObsPrimary","page":"Obs","title":"ObsPrimary","text":"","category":"section"},{"location":"obs/index_obs.html","page":"Obs","title":"Obs","text":" Corr\n corr_obs\n comp_t0\n comp_fvc","category":"page"},{"location":"obs/index_obs.html#HVPobs.Obs.Corr","page":"Obs","title":"HVPobs.Obs.Corr","text":"Corr(a::Vector{uwreal}, id::String, gamma::String)\n\nThe struct Corr stores data, ensemble id and gamma structure for two and three-point correlation functions. \n\n\n\n\n\n","category":"type"},{"location":"obs/index_obs.html#HVPobs.Obs.corr_obs","page":"Obs","title":"HVPobs.Obs.corr_obs","text":"corr_obs(cd::CData; real::Bool=true, rw::Union{Array{Float64,2}, Vector{Array{Float64,2}}, Nothing}=nothing, L::Int64=1, nms::Union{Int64, Nothing}=nothing)\n\nCreates a Corr struct with the given CData struct cdata  (read from read_hvp_data or read_mesons_data) for a single or multiple replicas. The flag real select the real or imaginary part of the correlator. If rw is specified, the method applies reweighting. rw is passed as a matrix of Float64 (read_ms1) The correlator can be normalized with the volume factor if L is fixed. The flag nms can be used to specify the maximum number of existing configurations in a given ensemble.\n\n\ndata = read_hvp_data(path, id)\nrw = read_ms1(path_rw)\n\ncorr = corr_obs(data)\ncorr_r = corr_obs(data, rw=rw)\n\n\n\n\n\n","category":"function"},{"location":"obs/index_obs.html#HVPobs.Obs.comp_t0","page":"Obs","title":"HVPobs.Obs.comp_t0","text":"comp_t0(Y::YData, plat::Vector{Int64}; L::Int64, pl::Bool=false, rw::Union{Matrix{Float64}, Nothing}=nothing, npol::Int64=2, wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing, info::Bool=false)\n\ncomp_t0(Y::Vector{YData}, plat::Vector{Int64}; L::Int64, pl::Bool=false, rw::Union{Vector{Matrix{Float64}}, Nothing}=nothing, npol::Int64=2, wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing, info::Bool=false)\n\nComputes t0 using the energy density of the action Ysl(Yang-Mills action). t0 is computed in the plateau plat. A polynomial interpolation in t is performed to find t0, where npol is the degree of the polynomial (linear fit by default)\n\nThe flag pl allows to show the plot. \n\nThe flag info provides extra output that contains information about the primary observables. The function returns the primary observables WY and W (it returns the observable <Y> if rw=nothing)\n\n#Single replica\nY = read_ms(path)\nrw = read_ms(path_rw)\n\nt0, Yobs = comp_t0(Y, [38, 58], L=32, info=true)\nt0_r, WYobs, Wobs = comp_t0(Y, [38, 58], L=32, rw=rw, info=true)\n\n#Two replicas\nY1 = read_ms(path1)\nY2 = read_ms(path2)\nrw1 = read_ms(path_rw1)\nrw2 = read_ms(path_rw2)\n\nt0 = comp_t0([Y1, Y2], [38, 58], L=32, pl=true)\nt0_r = comp_t0(Y, [38, 58], L=32, rw=[rw1, rw2], pl=true)\n\n\n\n\n\n\n","category":"function"},{"location":"obs/index_obs.html#HVPobs.Obs.comp_fvc","page":"Obs","title":"HVPobs.Obs.comp_fvc","text":"comp_fvc(path::String)\n\nComputes the finite volume corrections from the jackknife samples and returns the uwreal object as an Array{uwreal}(n2max, Thalf),  where n2max is the number (squared) of pion wrapping around the torus and Thalf is half of the lattice temporal extent.  The data file to be passed are typically named fvcorrblatgslrng0000err3ensIDwin0.dat \n\npath = /path/to/fv_corr_blat_gslrng0000_err3_ensID.dat\ndata = comp_fvc(path)\n\n\n\n\n\n","category":"function"},{"location":"obs/index_obs.html#ObsSpectrum","page":"Obs","title":"ObsSpectrum","text":"","category":"section"},{"location":"obs/index_obs.html","page":"Obs","title":"Obs","text":"meff","category":"page"},{"location":"obs/index_obs.html#HVPobs.Obs.meff","page":"Obs","title":"HVPobs.Obs.meff","text":"meff(obs::Vector{uwreal}, plat::Vector{Int64}; pl::Bool=false, data::Bool=false, wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)\nmeff(corr::Corr, plat::Vector{Int64}; pl::Bool=true, data::Bool=false, wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)\n\nComputes effective mass for a given correlator corr at a given plateau plat. Correlator can be passed as an Corr struct or Vector{uwreal}.\n\nThe flags pl and data allow to show the plots and return data as an extra result.\n\ndata = read_mesons(path, \"G5\", \"G5\")\ncorr_pp = corr_obs.(data)\nm = meff(corr_pp[1], [50, 60], pl=false)\n\n\n\n\n\n","category":"function"},{"location":"obs/index_obs.html#ObsTools","page":"Obs","title":"ObsTools","text":"","category":"section"},{"location":"obs/index_obs.html","page":"Obs","title":"Obs","text":"plat_av\nfrwd_bckwrd_symm!","category":"page"},{"location":"obs/index_obs.html#HVPobs.Obs.plat_av","page":"Obs","title":"HVPobs.Obs.plat_av","text":"@doc raw     plat_av(obs::Vector{uwreal}, plat::Vector{Int64}; wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)\n\nThis function performs the plateau average of obs over the plateau plat. The wpm flag can be used to control error analysis according to the ADerrors documentation.\n\n\n\n\n\n","category":"function"},{"location":"obs/index_obs.html#HVPobs.Obs.frwd_bckwrd_symm!","page":"Obs","title":"HVPobs.Obs.frwd_bckwrd_symm!","text":"@doc raw     frwdbckwrdsymm!(obs::Vector{uwreal})     frwdbckwrdsymm!(corr::Corr)\n\nThis function performs the forward-backward symmetrization of a given obs passed as input.\n\n\n\n\n\n","category":"function"},{"location":"obs/index_obs.html#ObsImprovement","page":"Obs","title":"ObsImprovement","text":"","category":"section"},{"location":"obs/index_obs.html","page":"Obs","title":"Obs","text":"improve_corr_vkvk!\nimprove_corr_vkvk_cons!\nZV\nZV_set2\ncv_loc\ncv_loc_set2\ncv_cons\ncv_cons_set2\nbv\nbv_set2\nbv_bar\nbv_bar_set2","category":"page"},{"location":"obs/index_obs.html#HVPobs.Obs.improve_corr_vkvk!","page":"Obs","title":"HVPobs.Obs.improve_corr_vkvk!","text":"improve_corr_vkvk!(vkvk::Vector{uwreal}, vkt0k::Vector{uwreal}, cv::Union{Float64,uwreal}; std::Bool=false)\nimprove_corr_vkvk!(vkvk::Corr, t0tk::Corr, cv::Union{Float64,uwreal}; std::Bool=false)\n\nGiven the vector-vector vkvkand vector-tensor vkt0k correlators as input, together with the cv improvement coefficient, this function overwrite the vkvk correlator  with the Symanzik improved version.  The flag std controls the derivative discretization of the vkt0k correlators. If false, it uses the improved version of the derivative, if true it uses the standard symmetric derivative.\n\n```@example\ncdata_vv = read_hvp_data(path_to_vv, id)\ncdata_vt = read_hvp_data(path_to_vt, id)\n\nrw = read_ms1(path, v=\"1.4\")\n\ncorr_vv = corr_obs(cdata_vv, rw=rw)\ncorr_vt = corr_obs(cdata_vt, rw=rw)\n\ncv = cv_loc(beta)\nimprove_corr_vkvk!(corr_vv, corr_vt, cv )\n```\n\n\n\n\n\n","category":"function"},{"location":"obs/index_obs.html#HVPobs.Obs.improve_corr_vkvk_cons!","page":"Obs","title":"HVPobs.Obs.improve_corr_vkvk_cons!","text":"improve_corr_vkvk_cons!(vkvk::Vector{uwreal}, vkt0k_l::Vector{uwreal}, vkt0k_c::Vector{uwreal}, cv_l::Union{Float64,uwreal}, cv_c::Union{Float64,uwreal}; std::Bool=false)\nimprove_corr_vkvk_cons!(vkvk::Corr, t0tk_l::Corr, t0tk_c::Corr, cv_l::Union{Float64,uwreal}, cv_c::Union{Float64,uwreal}; std::Bool=false)\n\nGiven the local-conserved vector-vector vkvk, the local-conserved vector-tensor vkt0k and the local-local vector-tensor vkt0k correlators as input, together with the cvlocal and cvcons improvement coefficient, this function overwrite the local-conserved vkvk correlator  with the Symanzik improved version.  The flag std controls the derivative discretization of the vkt0k correlators. If false, it uses the improved version of the derivative, if true it uses the standard symmetric derivative.\n\n```@example\ncdata_vv = read_hvp_data(path_to_vv, id)\ncdata_vt_loc = read_hvp_data(path_to_vt_loc, id)\ncdata_vt_cons = read_hvp_data(path_to_vt_cons, id)\n\nrw = read_ms1(path, v=\"1.4\")\n\ncorr_vv = corr_obs(cdata_vv, rw=rw)\ncorr_vt_loc = corr_obs(cdata_vt_loc, rw=rw)\ncorr_vt_cons = corr_obs(cdata_vt_cons, rw=rw)\n\ncv_loc  = cv_loc(beta)\ncv_cons = cv_cons(beta)\n\nimprove_corr_vkvk!(corr_vv, corr_vt_loc, corr_vt_cons, cv_loc, cv_cons)\n```\n\n\n\n\n\n","category":"function"},{"location":"obs/index_obs.html#HVPobs.Obs.ZV","page":"Obs","title":"HVPobs.Obs.ZV","text":"ZV(beta::Float64)\n\nGiven the coupling beta, this function returns the ZV renormalization constant for the vector current based on the results of 1811.08209. \n\n\n\n\n\n","category":"function"},{"location":"obs/index_obs.html#HVPobs.Obs.ZV_set2","page":"Obs","title":"HVPobs.Obs.ZV_set2","text":"ZV_set2(beta::Float64)\n\nGiven the coupling beta, this function returns the ZV renormalization constant for the vector current based on the SF results of 2010.09539. \n\n\n\n\n\n","category":"function"},{"location":"obs/index_obs.html#HVPobs.Obs.cv_loc","page":"Obs","title":"HVPobs.Obs.cv_loc","text":"cv_loc(beta::Float64)\n\nGiven the coupling beta, this function returns the cv improvement coefficient  for the local vector current based on the results of  1811.08209. \n\n\n\n\n\n","category":"function"},{"location":"obs/index_obs.html#HVPobs.Obs.cv_loc_set2","page":"Obs","title":"HVPobs.Obs.cv_loc_set2","text":"cvlocset2(beta::Float64)\n\nGiven the coupling beta, this function returns the cv improvement coefficient  for the local vector current based on the SF results of 2010.09539. \n\n\n\n\n\n","category":"function"},{"location":"obs/index_obs.html#HVPobs.Obs.cv_cons","page":"Obs","title":"HVPobs.Obs.cv_cons","text":"cv_cons(beta::Float64)\n\nGiven the coupling beta, this function returns the cv improvement coefficient  for the conserved vector current based on the results of  1811.08209. \n\n\n\n\n\n","category":"function"},{"location":"obs/index_obs.html#HVPobs.Obs.cv_cons_set2","page":"Obs","title":"HVPobs.Obs.cv_cons_set2","text":"cvconsset2(beta::Float64)\n\nGiven the coupling beta, this function returns the cv improvement coefficient  for the conserved vector current based on the SF results of 2010.09539. \n\n\n\n\n\n","category":"function"},{"location":"obs/index_obs.html#HVPobs.Obs.bv","page":"Obs","title":"HVPobs.Obs.bv","text":"bv(beta::Float64)\n\nGiven the coupling beta, this function returns the bv improvement coefficient for  the vector carrent based on the results of 1811.08209. \n\n\n\n\n\n","category":"function"},{"location":"obs/index_obs.html#HVPobs.Obs.bv_set2","page":"Obs","title":"HVPobs.Obs.bv_set2","text":"bv_set2(beta::Float64)\n\nGiven the coupling beta, this function returns the bv coefficient for the  improvement of the ZV renormalization constants. Values from 1805.07401\n\n\n\n\n\n","category":"function"},{"location":"obs/index_obs.html#HVPobs.Obs.bv_bar","page":"Obs","title":"HVPobs.Obs.bv_bar","text":"bv_bar(beta::Float64)\n\nGiven the coupling beta, this function returns the bv_bar improvement coefficient for  the vector carrent based on the results of 1811.08209. \n\n\n\n\n\n","category":"function"},{"location":"obs/index_obs.html#HVPobs.Obs.bv_bar_set2","page":"Obs","title":"HVPobs.Obs.bv_bar_set2","text":"bvbarset2(beta::Float64)\n\nGiven the coupling beta, this function returns the bv_bar coefficient for the  improvement of the Zv renormalization constants. Values from 1805.07401\n\n\n\n\n\n","category":"function"},{"location":"api/index_api.html","page":"API","title":"API","text":"","category":"page"},{"location":"index.html#HVPobs-Documentation","page":"Home","title":"HVPobs Documentation","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"A Julia tool for the analysis of LatticeQCD  Hadronic Vacuum Polarizations for the g-2. ","category":"page"},{"location":"index.html#Contents","page":"Home","title":"Contents","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Pages = [\n        \"data/index_data.md\",\n        \"obs/index_obs.md\", \n        \"fits/index_fits.md\"\n        ]\nDepth = 3","category":"page"}]
}