@doc raw"""
    improve_corr_vkvk!(vkvk::Vector{uwreal}, vkt0k::Vector{uwreal}, cv::Union{Float64,uwreal}; std::Bool=false)
    improve_corr_vkvk!(vkvk::Corr, t0tk::Corr, cv::Union{Float64,uwreal}; std::Bool=false)

Given the vector-vector `vkvk`and vector-tensor `vkt0k` correlators as input, together with the cv improvement coefficient,
this function overwrite the `vkvk` correlator  with the Symanzik improved version. 
The flag `std` controls the derivative discretization of the `vkt0k` correlators. If false, it uses the improved version of the derivative,
if true it uses the standard symmetric derivative.

    ```@example
    cdata_vv = read_hvp_data(path_to_vv, id)
    cdata_vt = read_hvp_data(path_to_vt, id)

    rw = read_ms1(path, v="1.4")

    corr_vv = corr_obs(cdata_vv, rw=rw)
    corr_vt = corr_obs(cdata_vt, rw=rw)

    cv = cv_loc(beta)
    improve_corr_vkvk!(corr_vv, corr_vt, cv )
    ```
"""
function improve_corr_vkvk!(vkvk::Vector{uwreal}, vkt0k::Vector{uwreal}, cv::Union{Float64,uwreal}; std::Bool=false, treelevel::Bool=false)
    
    der_t0tk = improve_derivative(vkt0k, std=std, treelevel=treelevel)
    vkvk[2:end] = vkvk[2:end] .+ cv .* der_t0tk
end
improve_corr_vkvk!(vkvk::Corr, t0tk::Corr, cv::Union{Float64,uwreal}; std::Bool=false, treelevel::Bool=false) = improve_corr_vkvk!(vkvk.obs, t0tk.obs, cv, std=std, treelevel=treelevel)

@doc raw"""
    improve_corr_vkvk_cons!(vkvk::Vector{uwreal}, vkt0k_l::Vector{uwreal}, vkt0k_c::Vector{uwreal}, cv_l::Union{Float64,uwreal}, cv_c::Union{Float64,uwreal}; std::Bool=false)
    improve_corr_vkvk_cons!(vkvk::Corr, t0tk_l::Corr, t0tk_c::Corr, cv_l::Union{Float64,uwreal}, cv_c::Union{Float64,uwreal}; std::Bool=false)

Given the local-conserved vector-vector `vkvk`, the local-conserved vector-tensor `vkt0k` and the local-local vector-tensor `vkt0k` correlators as input, together with the cv_local and cv_cons improvement coefficient,
this function overwrite the local-conserved `vkvk` correlator  with the Symanzik improved version. 
The flag `std` controls the derivative discretization of the `vkt0k` correlators. If false, it uses the improved version of the derivative,
if true it uses the standard symmetric derivative.

    ```@example
    cdata_vv = read_hvp_data(path_to_vv, id)
    cdata_vt_loc = read_hvp_data(path_to_vt_loc, id)
    cdata_vt_cons = read_hvp_data(path_to_vt_cons, id)

    rw = read_ms1(path, v="1.4")

    corr_vv = corr_obs(cdata_vv, rw=rw)
    corr_vt_loc = corr_obs(cdata_vt_loc, rw=rw)
    corr_vt_cons = corr_obs(cdata_vt_cons, rw=rw)

    cv_loc  = cv_loc(beta)
    cv_cons = cv_cons(beta)
    
    improve_corr_vkvk!(corr_vv, corr_vt_loc, corr_vt_cons, cv_loc, cv_cons)
    ```
"""
function improve_corr_vkvk_cons!(vkvk::Vector{uwreal}, vkt0k_l::Vector{uwreal}, vkt0k_c::Vector{uwreal}, cv_l::Union{Float64,uwreal}, cv_c::Union{Float64,uwreal}; std::Bool=false, treelevel::Bool=false)
    der_t0tk_l = improve_derivative(vkt0k_l, std=std, treelevel=treelevel)
    der_t0tk_c = improve_derivative(vkt0k_c, std=std, treelevel=treelevel)

    vkvk[2:end] = vkvk[2:end] .+ cv_l .* der_t0tk_c .+ cv_c .* der_t0tk_l
end
improve_corr_vkvk_cons!(vkvk::Corr, t0tk_l::Corr, t0tk_c::Corr, cv_l::Union{Float64,uwreal}, cv_c::Union{Float64,uwreal}; std::Bool=false, treelevel::Bool=false) = improve_corr_vkvk_cons!(vkvk.obs, t0tk_l.obs, t0tk_c.obs, cv_l, cv_c, std=std, treelevel=treelevel)

function improve_derivative(corr::Vector{uwreal}; std::Bool=false, treelevel::Bool=false)
    if std
        dcorr = (corr[3:end] - corr[1:end-2]) / 2 
        dcorr[1] = corr[3] - corr[2]
        push!(dcorr, corr[end] - corr[end-1])
        return dcorr
    else
        if treelevel
            tvals = Float64.(collect(0:length(corr)-1))
        else
            tvals = vcat(collect(0:length(corr)/2), -reverse(collect(1:length(corr)/2-1))...)
        end
        corr_aux = corr .* tvals.^2
        dcorr_aux = improve_derivative(corr_aux, std=true)
        dcorr = 1 ./ tvals[2:end].^2 .* (dcorr_aux .- 2 .* tvals[2:end] .* corr[2:end])
        return dcorr
    end
end


@doc raw"""
Get the perturbatively estimated value for the improvement coefficients cV

"""
function cv_pert(beta::Float64)
    g0sq = 6 ./beta
    CF = 4 ./3
    cV1=-0.01030*4/3
    return cV1 * g0sq
end
@doc raw"""
ZV(beta::Float64)

Given the coupling beta, this function returns the ZV renormalization constant for
the vector current based on the results of 1811.08209. 
"""
function ZV(beta::Float64)
    g2 = 6/beta
    # CC = [[1.11887e-8, -1.34335e-4*1.11887e-4*1.61288e-4, 4.02275e-5*1.11887e-4*1.44633e-5] [-1.34335e-4*1.11887e-4*1.61288e-4, 1.61288e-8, -4.82986e-5*1.44633e-5*1.61288e-4] [4.02275e-5*1.11887e-4*1.44633e-5, -4.82986e-5*1.44633e-5*1.61288e-4, 1.44633e-10 ]]
    CC = [[1.31619, 4.92750 , 6.15758] [4.92750, 66.8321, 75.3218] [ 6.15758,  75.3218, 85.2733 ]] .*1e-6
    p = cobs([-0.2542, -0.0961, -0.4796], CC, [1,2,3])
    zv = 1 - 0.10057*g2 *(1 + p[1]*g2 + p[2]*g2^2)/(1 +p[3]*g2)
    return zv
end

@doc raw"""
bv(beta::Float64)

Given the coupling beta, this function returns the bv improvement coefficient for 
the vector carrent based on the results of 1811.08209. 
"""
function bv(beta::Float64)

    g2 = 6/beta
    CC = [[36.7139, 12.6698] [12.6698, 4.41224]] .*1e-4
    p = cobs([-0.184, -0.444], CC, [4,5] )
    bvv = 1 + 0.11813 * g2 * (1 + p[1]*g2) / ( 1 + p[2]*g2)
    return bvv
end

@doc raw"""
bv_bar(beta::Float64)

Given the coupling beta, this function returns the bv_bar improvement coefficient for 
the vector carrent based on the results of 1811.08209. 
"""
function  bv_bar(beta::Float64)

    g2 = 6/beta
    CC = [[1.061463, 14.53004 ] [14.53004, 248.5266]] .*1e-8
    p = cobs([0.00112, -0.5577], CC, [6,7] )
    bvbar = ( p[1]*g2^2) / ( 1 + p[2]*g2)
    return bvbar
end

@doc raw"""
cv_loc(beta::Float64)

Given the coupling beta, this function returns the cv improvement coefficient 
for the local vector current based on the results of  1811.08209. (updated values not yet published)
"""
function cv_loc(beta::Float64)
    g2 = 6 / beta
    b0 = 9 / (4*π)^2
    p = [-2603.598964596533, 2097.241663044384] 
    return -0.01030 * 4 / 3 * g2 * (1 + exp(-1/(2*b0*g2)) * (p[1] + p[2]*g2) )
end
function cv_loc_old(beta::Float64)
    # Old values of Mainz improvement coefficients
    g2 = 6/beta
    p = uwreal([0.15, 0.35], "cv_loc 1-loop")
    cv = -0.01030 * 4 / 3 * g2 * (1. + value(p)*g2)
    return cv
end

@doc raw"""
cv_cons(beta::Float64)

Given the coupling beta, this function returns the cv improvement coefficient 
for the conserved vector current based on the results of  1811.08209. (updated values not yet published)
"""
function cv_cons(beta::Float64)
    g2 = 6 / beta
    b0 = 9 / (4*π)^2
    p = [0.024348066540864678, 321.0689682507245] 
    return 0.5 - p[1] * g2 * (1 + exp(-1/(2*b0*g2)) * (p[2]*g2) )
end
function cv_cons_old(beta::Float64)
   # Old values of Mainz improvement coefficients
    g2 = 6/beta
    p = uwreal([0.093, 0.13], "cv_cons 1-loop")
    cv = 0.5 * (1 - value(p)*g2)
    return cv
end

# SET 2 of imrovement coefficients
@doc raw"""
cv_loc_set2(beta::Float64)

Given the coupling beta, this function returns the cv improvement coefficient 
for the local vector current based on the SF results of 2010.09539. 
"""
function cv_loc_set2(beta::Float64)
    # g2 = 6/beta
    # CC = [[1.97358e6^2, -1.15750e6*1.97358e6*6.78869e5 ] [-1.15750e6*1.97358e6*6.78869e5, 6.78869e5^2]]
    # p = cobs([-3.039e3, 2.649e3], CC, [8,9])
    # cv = -0.01030*4/3*g2 * (1+ exp(-1/(2*(9/(16*pi^2))*g2))*(-3.039e3 +2.649e3*g2))
    cv = Dict(
        3.34 => -0.346,
        3.4  => -0.299,
        3.46 => -0.259,
        3.55 => -0.209,
        3.7  => -0.147, 
        3.85  => -0.105 
    )
    return cv[beta]
end

@doc raw"""
cv_cons_set2(beta::Float64)

Given the coupling beta, this function returns the cv improvement coefficient 
for the conserved vector current based on the SF results of 2010.09539. 
"""
function cv_cons_set2(beta::Float64)
    # g2 = 6/beta
    # cv = 0.5 - 2.112e-2 * g2 * (1 + exp(-1/(2*(9/(16*pi^2))*g2)) * 506.5*g2 )
    cv = Dict(
        3.34 => 0.201,
        3.4  => 0.232,
        3.46 => 0.259,
        3.55 => 0.294,
        3.7  => 0.340, 
        3.85 => 0.374 
    )
    return cv[beta]
end

@doc raw"""
ZV_set2(beta::Float64)

Given the coupling beta, this function returns the ZV renormalization constant for
the vector current based on the SF results of 2010.09539. 
"""
function ZV_set2(beta::Float64)
    g2 = 6/beta
    CC = [[1.11887e-8, -1.34335e-8, 4.02275e-10] [-1.34335e-8, 1.61288e-8, -4.82986e-10] [4.02275e-10, -4.82986e-10, 1.44633e-10 ]]
    p = cobs([-1.907e-1, 2.199e-1, -7.492e-2], CC, [8,9,10])
    zv = 1 - 0.075427*4/3*g2 + p[1]*g2^2 + p[2]*g2^3 +p[3]*g2^4
    return zv
end

@doc raw"""
bv_set2(beta::Float64)

Given the coupling beta, this function returns the bv coefficient for the 
improvement of the ZV renormalization constants.
Values from 1805.07401
"""
function bv_set2(beta::Float64)
    g2 = 6/beta
    CC = [[86.1424, -20.5770, 23.7151, 78.9559] [-20.5770, 8.1121, -6.9724, -17.7160] [23.7151, -6.9724, 7.1199, 21.3673] [78.9559, -17.7160, 21.3673, 72.9602]] .* 1e-5
    p = cobs([-0.40040, 0.04352, -0.03852, -0.48803], CC, [11,12,13,14])
    bv = (1+ p[1]*g2 + p[2]*g2^2 + p[3]*g2^3) / (1 + p[4]*g2)
    return bv
end

@doc raw"""
bv_bar_set2(beta::Float64)

Given the coupling beta, this function returns the bv_bar coefficient for the 
improvement of the Zv renormalization constants.
Values from 1805.07401
"""
function bv_bar_set2(beta::Float64)
    g2 = 6/beta
    CC = [[4.9676, -4.8982, 3.0473, 3.2513] [-4.8982, 6.4514, -3.6515, -2.6826] [3.0473, -3.6515, 2.1830, 1.8785] [3.2513, -2.6826, 1.8785, 2.4558]] .*1e-5
    p = cobs([-0.43101, 0.04109, -0.03911, -0.51771], CC, [15,16,17,18])
    rvbar = (1 + p[1]*g2 + p[2]*g2^2 + p[3]*g2^3) / (1 + p[4]*g2)
    bvbar = (rvbar - bv_set2(beta)) / 3
    return bvbar

end
@doc raw"""
    ca(beta::Float64)
Given the coupling beta, this function returns the ca coefficient for the 
improvement of the axial correlator. Taken from 1502.04999
"""
function ca(beta::Float64)
    g2 = 6/beta
    ca = -0.006033 * g2 * ( 1 + exp(9.2056 - 13.9847/g2))
    return ca
end

@doc raw"""
Za_l_sub(beta::Float64)
Given the coupling beta, this function returns the renormalization factor Za^l sub from LCP-2 of 1808.09236  
This formula is valid up to β=3.85
"""
function Za_l_sub(beta::Float64)
    g2 = 6/beta
    CC = [[0.229571866, -0.278151898, 0.084105454] [-0.278151898, 0.337131945, -0.101975449] [0.084105454, -0.101975449, 0.030856380]] .* 1e-1
    p = cobs([1.35510, -0.501106, 0.091656], CC, [19,20,21])
    Za = p[1] + p[2]*g2 + p[3]*g2^2
    return Za 
end

@doc raw"""
ZP(beta::Float64)
Given the coupling beta, this function returns the renormalization factor ZP from 1802.05243.
"""
function ZP(beta::Float64)

    C = [[0.375369e-6, 0.429197e-6, -0.186896e-5] [0.429197e-6, 0.268393e-4, 0.686776e-4] [-0.186896e-5, 0.686776e-4, 0.212386e-3]]
    z = cobs([0.348629, 0.020921, 0.070613], C, [22,23,24])
    zp = z[1] + z[2] * (beta - 3.79) + z[3] * (beta - 3.79)^2
    return zp
end

@doc raw"""
    ba_minus_bp(beta::Float64)
Given the coupling beta, this function returns the ba-bp coefficient for the 
improvement of the ZA/ZP renormalization constants.
Values from  hep-lat/9806015
"""
function ba_minus_bp(beta::Float64)
    g2 = 6/beta
    # from scale setting Bruno et al. 
    res = - 0.0012*g2
    return res
end
@doc raw"""
    ba_imp(beta::Float64)
Given the coupling beta, this function returns the ba coefficient for the 
improvement of the Za renormalization constants
Values from one-loop computation hep-lat/9806015
"""
function ba_imp(beta::Float64)
    g2 = 6/beta
    ba = 1 + 0.1141*4/3  *g2
    return ba
end
