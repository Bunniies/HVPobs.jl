@doc raw"""
    corrConnected(path_data::String, ens::EnsInfo, sector::String; path_rw=Union{Nothing,String}=nothing, impr::Bool=true, impr_set::String="1", std::Bool=false, L::Int64=1)

Computes the Vector-Vector  local-local and local-conserved connected correlators for a given EnsInfo `ens` and sector `sector`.
The supported sectors are: light, strange, charm, charm_plus.
It returns the forward-backward symmetrised local-local and local-conserved correlators as Corr objects.
No charge factor is included in the computation. 
A sign flip is performed on the correlator according to Eq 16 of 2206.06582.


Optional flags:    

    - path_rw  : if !isnothing(path_rw), correlators are reweighted.  
    - impr     : if true, the correlator is improved.  
    - impr_set : either "1" or "2", to select set1 (Mainz) or set2 (ALPHA) improvement coefficients accordingly.  
    - std      : if true, standard symmetric derivatives are used for the vector-tensor correlator. If false, improved derivatives are used.
    - L        : correlators are normalised with the volume L^3. L=1 by default.    

Examples:
```@example
ens = EnsInfo("H101")
gvv_ll, gvv_lc = corrConnected(pathToData, ens, "ligth", path_rw=pathToRwf, impr=true, impr_set="1", std=false, L=1) 
```
"""
function corrConnected(path_data::String, ens::EnsInfo, sector::String; path_rw=Union{Nothing,String}=nothing, impr::Bool=true, impr_set::String="1", std::Bool=false, L::Int64=1)
    if sector ∉ ["light", "strange", "charm", "charm_plus"]
        error("Flavour sector '$(sector)' not recognised.")
    end
    # local-local VV
    Gamma_l = ["V1V1", "V2V2", "V3V3", "V1T10", "V2T20", "V3T30"]
    
    v1v1 = get_corr(path_data, ens, sector, Gamma_l[1], path_rw=path_rw, frw_bcwd=false, L=L)
    v2v2 = get_corr(path_data, ens, sector, Gamma_l[2], path_rw=path_rw, frw_bcwd=false, L=L)
    v3v3 = get_corr(path_data, ens, sector, Gamma_l[3], path_rw=path_rw, frw_bcwd=false, L=L)
    gvv_ll = Corr(-1 .* (v1v1.obs .+ v2v2.obs .+ v3v3.obs)./3, ens.id, "Gvv_ll_"*sector)
        
    # local-conserved VV
    Gamma_c = ["V1V1c", "V2V2c", "V3V3c", "V1cT10", "V2cT20", "V3cT30"]
    v1v1_c = get_corr(path_data, ens, sector, Gamma_c[1], path_rw=path_rw, frw_bcwd=false, L=L)
    v2v2_c = get_corr(path_data, ens, sector, Gamma_c[2], path_rw=path_rw, frw_bcwd=false, L=L)
    v3v3_c = get_corr(path_data, ens, sector, Gamma_c[3], path_rw=path_rw, frw_bcwd=false, L=L)
    gvv_lc = Corr(-1 .* (v1v1_c.obs .+ v2v2_c.obs .+ v3v3_c.obs)./3, ens.id, "Gvv_lc_"*sector)

    
    if impr
        # local-local VT
        v1t10 = get_corr(path_data, ens, sector, Gamma_l[4], path_rw=path_rw, frw_bcwd=false, L=L)
        v2t20 = get_corr(path_data, ens, sector, Gamma_l[5], path_rw=path_rw, frw_bcwd=false, L=L)
        v3t30 = get_corr(path_data, ens, sector, Gamma_l[6], path_rw=path_rw, frw_bcwd=false, L=L)
        gvt_ll = Corr(-1 .* (v1t10.obs .+ v2t20.obs .+ v3t30.obs)./3, ens.id, "Gvt_ll_"*sector)

        # local-consereved VT
        v1t10_c = get_corr(path_data, ens, sector, Gamma_c[4], path_rw=path_rw, frw_bcwd=false, L=L)
        v2t20_c = get_corr(path_data, ens, sector, Gamma_c[5], path_rw=path_rw, frw_bcwd=false, L=L)    
        v3t30_c = get_corr(path_data, ens, sector, Gamma_c[6], path_rw=path_rw, frw_bcwd=false, L=L)
        gvt_lc = Corr(-1 .* (v1t10_c.obs .+ v2t20_c.obs .+ v3t30_c.obs)./3, ens.id, "Gvt_lc_"*sector)

        beta = ens.beta
        if impr_set == "1"
            cv_l = cv_loc(beta)
            cv_c = cv_cons(beta)
        elseif impr_set == "1old"
            cv_l = cv_loc_old(beta)
            cv_c = cv_cons_old(beta)
        elseif impr_set =="2"
            cv_l = cv_loc_set2(beta)
            cv_c = cv_cons_set2(beta)
        end

        improve_corr_vkvk!(gvv_ll, gvt_ll, 2*cv_l, std=std)
        improve_corr_vkvk_cons!(gvv_lc, gvt_ll, gvt_lc, cv_l, cv_c, std=std)
    end
    frwd_bckwrd_symm!(gvv_ll)
    frwd_bckwrd_symm!(gvv_lc)

    return gvv_ll, gvv_lc
end

@doc raw"""
    corrDisconnected(path_data::String, ens::EnsInfo, fl::String; path_rw::Union{Nothing,String}=nothing, impr::Bool=true, impr_set::String="1", std::Bool=true, L::Int64=1)

Computes the Vector-Vector  local-local, local-conserved and consereved-conserved disconnected correlators for a given EnsInfo `ens` and flavour `fl`.
The supported flavours are: [08, 0c, 80, 88, 8c, c0, c8, cc].
It returns the forward-backward symmetrised local-local, local-conserved, conserved-conserved correlators as Corr objects.
No charge factor is included in the computation. 
A sign flip is performed on the correlator according to Eq 16 of 2206.06582.


Optional flags:    

    - path_rw  : if !isnothing(path_rw), correlators are reweighted.  
    - impr     : if true, the correlator is improved.  
    - impr_set : either "1" or "2", to select set1 (Mainz) or set2 (ALPHA) improvement coefficients accordingly.  
    - std      : if true, standard symmetric derivatives are used for the vector-tensor correlator. If false, improved derivatives are used.
    - L        : correlators are normalised with the volume L^3. L=1 by default.    

Examples:
```@example
ens = EnsInfo("N200")
gvv_ll, gvv_lc , gvv_cc = corrDonnected(pathToData, ens, "cc", path_rw=pathToRwf, impr=true, impr_set="1", std=false, L=1) 
```
"""
function corrDisconnected(path_data::String, ens::EnsInfo, fl::String; path_rw::Union{Nothing,String}=nothing, impr::Bool=true, impr_set::String="1", std::Bool=true, L::Int64=1)
    if fl ∉ ["08", "0c", "80", "88", "8c", "c0", "c8", "cc"]
        error("Unrecognised flavour structure $(fl): choose from [08, 0c, 80, 88, 8c, c0, c8, cc]")
    end

    if fl ∈ ["08", "80"]
        return corrDisconnected80(path_data, ens, path_rw=path_rw, impr=impr, impr_set=impr_set, std=std)
    end

    corr_dict = get_corr_disc(path_data, ens, fl, path_rw=path_rw, frw_bcwd=false, L=L)
    g_ll = Corr(corr_dict["VV"].obs, ens.id,   "G"*fl*"_ll_disc")
    g_lc = Corr(corr_dict["VVc"].obs, ens.id,  "G"*fl*"_lc_disc")
    g_cc = Corr(corr_dict["VcVc"].obs, ens.id, "G"*fl*"_cc_disc")

    if impr
        beta = ens.beta
        if impr_set == "1"
            cv_l = cv_loc(beta) 
            cv_c = cv_cons(beta) 
        elseif impr_set == "1old"
            cv_l = cv_loc_old(beta)
            cv_c = cv_cons_old(beta)     
        elseif impr_set =="2"
            cv_l = cv_loc_set2(beta)
            cv_c = cv_cons_set2(beta)
        end
        improve_corr_vkvk!(g_ll, corr_dict["VT"], 2*cv_l, std=std)
        improve_corr_vkvk!(g_cc, corr_dict["VcT"], 2*cv_c, std=std)
        improve_corr_vkvk_cons!(g_lc, corr_dict["VT"], corr_dict["VcT"], cv_l, cv_c, std=std)
       
    end
   
    frwd_bckwrd_symm!(g_ll); renormalize!(g_ll, uwreal(-1.0))
    frwd_bckwrd_symm!(g_lc); renormalize!(g_lc, uwreal(-1.0))
    frwd_bckwrd_symm!(g_cc); renormalize!(g_cc, uwreal(-1.0))

    return g_ll, g_lc, g_cc
end

@doc raw"""
    corrDisconnected80(path_data::String, ens::EnsInfo; path_rw::Union{Nothing,String}=nothing, impr::Bool=true, impr_set::String="1", std::Bool=true, L::Int64=1)

Computes the Vector-Vector  local-local, local-conserved and consereved-conserved disconnected correlators for a given EnsInfo `ens` and flavour 08 or 80.
It returns the forward-backward symmetrised local-local, local-conserved, conserved-conserved correlators as Corr objects.
No charge factor is included in the computation. This function is called by corrDisconnected when flavour 08 or 80 are detected.
A sign flip is performed on the correlator according to Eq 16 of 2206.06582.

Optional flags:    

    - path_rw  : if !isnothing(path_rw), correlators are reweighted.  
    - impr     : if true, the correlator is improved.  
    - impr_set : either "1" or "2", to select set1 (Mainz) or set2 (ALPHA) improvement coefficients accordingly.  
    - std      : if true, standard symmetric derivatives are used for the vector-tensor correlator. If false, improved derivatives are used.
    - L        : correlators are normalised with the volume L^3. L=1 by default.    

Examples:
```@example
ens = EnsInfo("N200")
gvv_ll, gvv_lc , gvv_cc = corrDonnected80(pathToData, ens, "cc", path_rw=pathToRwf, impr=true, impr_set="1", std=false, L=1) 
```
"""
function corrDisconnected80(path_data::String, ens::EnsInfo; path_rw::Union{Nothing,String}=nothing, impr::Bool=true, impr_set::String="1", std::Bool=true, L::Int64=1)
    

    corr_dict_80 = get_corr_disc(path_data, ens, "80", path_rw=path_rw, frw_bcwd=false, L=L)
    corr_dict_08 = get_corr_disc(path_data, ens, "08", path_rw=path_rw, frw_bcwd=false, L=L)
    
    g_ll80 = Corr(corr_dict_80["VV"].obs, ens.id,   "G80_ll_disc")
    g_lc80 = Corr(corr_dict_80["VVc"].obs, ens.id,  "G80_lc_disc")
    g_cc80 = Corr(corr_dict_80["VcVc"].obs, ens.id, "G80_cc_disc")

    g_ll08 = Corr(corr_dict_08["VV"].obs, ens.id,   "G08_ll_disc")
    g_lc08 = Corr(corr_dict_08["VVc"].obs, ens.id,  "G08_lc_disc")
    g_cc08 = Corr(corr_dict_08["VcVc"].obs, ens.id, "G08_cc_disc")

    if impr
        beta = ens.beta
        if impr_set == "1"
            cv_l = cv_loc(beta) 
            cv_c = cv_cons(beta)      
        elseif impr_set == "1old"
            cv_l = cv_loc_old(beta)
            cv_c = cv_cons_old(beta)
        elseif impr_set =="2"
            cv_l = cv_loc_set2(beta)
            cv_c = cv_cons_set2(beta)
        end

        # local-local
        deriv_t8v0 = Obs.improve_derivative(corr_dict_80["TV"].obs, std=std) 
        deriv_v8t0 = Obs.improve_derivative(corr_dict_80["VT"].obs, std=std)

        deriv_t0v8 = Obs.improve_derivative(corr_dict_08["TV"].obs, std=std) 
        deriv_v0t8 = Obs.improve_derivative(corr_dict_08["VT"].obs, std=std)

        g_ll80.obs[2:end] = g_ll80.obs[2:end] + cv_l .* (-deriv_t8v0 + deriv_v8t0) 
        g_ll08.obs[2:end] = g_ll08.obs[2:end] + cv_l .* (-deriv_t0v8 + deriv_v0t8) 

        # local-conserved
        deriv_vc0t8 = Obs.improve_derivative(corr_dict_08["VcT"].obs, std=std)

        g_lc80.obs[2:end] = g_lc80.obs[2:end] + cv_l .* deriv_vc0t8 .+ cv_c .* deriv_v8t0

        # conserved-conserved
        deriv_t8vc0 = Obs.improve_derivative(corr_dict_80["TVc"].obs, std=std) 
        deriv_vc8t0 = Obs.improve_derivative(corr_dict_80["VcT"].obs, std=std)

        g_cc80.obs[2:end] = g_cc80.obs[2:end] + cv_c .* (-deriv_t8vc0 + deriv_vc8t0)
        
    end
    
    g_ll = Corr(g_ll80.obs , ens.id, "G08_ll_disc")
    g_lc = Corr(g_lc80.obs , ens.id, "G08_lc_disc")
    g_cc = Corr(g_cc80.obs , ens.id, "G08_cc_disc")

    frwd_bckwrd_symm!(g_ll); renormalize!(g_ll, uwreal(-1.0))
    frwd_bckwrd_symm!(g_lc); renormalize!(g_lc, uwreal(-1.0))
    frwd_bckwrd_symm!(g_cc); renormalize!(g_cc, uwreal(-1.0))

    return g_ll, g_lc, g_cc
end

@doc raw"""
    get_Z3(ens::EnsInfo; impr_set::String="1")

Given an EnsInfo `ens` and impr_set (Mainz -> 1, ALPHA -> 2), it returns the isovector G33 renormalization constant 
"""
function get_Z3(ens::EnsInfo; impr_set::String="1")
    
    ml = 0.5 * (1 / ens.kappa_l - 1 / ens.kappa_crit )
    ms = 0.5 * (1 / ens.kappa_s - 1 / ens.kappa_crit )
    trmq = 2*ml + ms
    if impr_set in ["1","1old"]
        Z3 = ZV(ens.beta) * (1. + bv_bar(ens.beta)*trmq + bv(ens.beta)*ml)
    elseif  impr_set == "2"
        Z3 = ZV_set2(ens.beta) * (1. + bv_bar_set2(ens.beta)*trmq + bv_set2(ens.beta)*ml)
    end
    return Z3
end

@doc raw"""
    get_Z8(ens::EnsInfo; impr_set::String="1")

Given an EnsInfo `ens` and impr_set (Mainz -> 1, ALPHA -> 2), it returns the isoscalar G88 renormalization constant 
"""
function get_Z8(ens::EnsInfo; impr_set::String="1")

    ml = 0.5 * (1 / ens.kappa_l - 1 / ens.kappa_crit )
    ms = 0.5 * (1 / ens.kappa_s - 1 / ens.kappa_crit )
    trmq = 2*ml + ms
    if impr_set in ["1","1old"]
        Z8 = ZV(ens.beta) * (1. + bv_bar(ens.beta)*trmq + bv(ens.beta)*(ml + 2*ms)/3)
    elseif impr_set == "2"
        Z8 = ZV_set2(ens.beta) * (1. + bv_bar_set2(ens.beta)*trmq + bv_set2(ens.beta)*(ml + 2*ms)/3)
    end
    return Z8
end

@doc raw"""
    get_Z08(ens::EnsInfo; impr_set::String="1")

Given an EnsInfo `ens` and impr_set (Mainz -> 1, ALPHA -> 2), it returns the  G08 renormalization constant 
"""
function get_Z08(ens::EnsInfo; impr_set::String="1")

    ml = 0.5 * (1 / ens.kappa_l - 1 / ens.kappa_crit )
    ms = 0.5 * (1 / ens.kappa_s - 1 / ens.kappa_crit )
    if impr_set in ["1","1old"]
        Z08 = 1/3 * ZV(ens.beta) * bv(ens.beta) * (2/sqrt(3)) * (ml - ms)
    elseif impr_set == "2"
        Z08 = 1/3 * ZV_set2(ens.beta) * bv_set2(ens.beta) * (2/sqrt(3)) * (ml - ms)
    end
    return Z08
end

@doc raw"""
    renormalize!(g::Vector{uwreal}, Z::uwreal)
    renormalize!(g::Corr, Z::uwreal)

Renormalise with Z a given Vector{uwreal} or  Corr object g. 
"""
function renormalize!(g::Vector{uwreal}, Z::uwreal)
    g[:] .*= Z
    return nothing
end
renormalize!(g::Corr, Z::uwreal) = renormalize!(g.obs, Z)