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
function improve_corr_vkvk!(vkvk::Vector{uwreal}, vkt0k::Vector{uwreal}, cv::Union{Float64,uwreal}; std::Bool=false)
    
    der_t0tk = improve_derivative(vkt0k, std=std)
    vkvk[2:end] = vkvk[2:end] .+ cv .* der_t0tk
end
improve_corr_vkvk!(vkvk::Corr, t0tk::Corr, cv::Union{Float64,uwreal}; std::Bool=false) = improve_corr_vkvk!(vkvk.obs, t0tk.obs, cv, std=std)

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
function improve_corr_vkvk_cons!(vkvk::Vector{uwreal}, vkt0k_l::Vector{uwreal}, vkt0k_c::Vector{uwreal}, cv_l::Union{Float64,uwreal}, cv_c::Union{Float64,uwreal}; std::Bool=false)
    
    der_t0tk_l = improve_derivative(vkt0k_l, std=std)
    der_t0tk_c = improve_derivative(vkt0k_c, std=std)

    vkvk[2:end] = vkvk[2:end] .+ cv_l .* der_t0tk_c .+ cv_c .* der_t0tk_l
end
improve_corr_vkvk_cons!(vkvk::Corr, t0tk_l::Corr, t0tk_c::Corr, cv_l::Union{Float64,uwreal}, cv_c::Union{Float64,uwreal}; std::Bool=false) = improve_corr_vkvk_cons!(vkvk.obs, t0tk_l.obs, t0tk_c.obs, cv_l, cv_c, std=std)

function improve_derivative(corr::Vector{uwreal}; std::Bool=false)
    if std
        dcorr = (corr[3:end] - corr[1:end-2]) / 2 
        dcorr[1] = corr[3] - corr[2]
        push!(dcorr, corr[end] - corr[end-1])
        return dcorr
    else
        tvals = Float64.(collect(0:length(corr)-1))
        corr_aux = corr .* tvals.^2
        dcorr_aux = (corr_aux[3:end] - corr_aux[1:end-2]) / 2
        dcorr_aux[1] = corr_aux[3] - corr_aux[2]
        push!(dcorr_aux, corr[end] - corr[end-1])
        dcorr = 1 ./ tvals[2:end].^2 .* (dcorr_aux .- 2 .* tvals[2:end] .* corr[2:end])
        return dcorr
    end
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
for the local vector current based on the results of  1811.08209. 
"""
function cv_loc(beta::Float64)
    
    g2 = 6/beta
    p = uwreal([0.15, 0.35], "cv_loc 1-loop")
    cv = -0.01030 * 4 / 3 * g2 * (1. + value(p)*g2)
    return cv
end

@doc raw"""
cv_cons(beta::Float64)

Given the coupling beta, this function returns the cv improvement coefficient 
for the conserved vector current based on the results of  1811.08209. 
"""
function cv_cons(beta::Float64)
   
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
    g2 = 6/beta
    #CC = [[1.97358e6^2, -1.15750e6*1.97358e6*6.78869e5 ] [-1.15750e6*1.97358e6*6.78869e5, 6.78869e5^2]]
    #p = cobs([-3.039e3, 2.649e3], CC, [8,9])
    cv = -0.01030*4/3*g2 * (1+ exp(-1/(2*(9/(16*pi^2))*g2))*(-3.039e3 +2.649e3*g2))
    return cv
end

@doc raw"""
cv_cons_set2(beta::Float64)

Given the coupling beta, this function returns the cv improvement coefficient 
for the conserved vector current based on the SF results of 2010.09539. 
"""
function cv_cons_set2(beta::Float64)
    g2 = 6/beta
    cv = 0.5 - 2.112e-2 * g2 * (1 + exp(-1/(2*(9/(16*pi^2))*g2)) * 506.5*g2 )
    return cv
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
    p = cobs([-0.43101, 0.04109, -0.03911, -0.51771], CC, [15,16,16,18])
    rvbar = (1 + p[1]*g2 + p[2]*g2^2 + p[3]*g2^3) / (1 + p[4]*g2)
    bvbar = (rvbar - bv_set2(beta)) / 3
    return bvbar

end