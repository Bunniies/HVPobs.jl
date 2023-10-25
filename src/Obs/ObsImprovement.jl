@doc raw"""
    improve_corr_vkvk!(vkvk::Vector{uwreal}, vkt0tk::Vector{uwreal}, cv::Float64)

Given vkvk and t0kt0k correlators as input, together with cv improvement coefficient,
this function improves the the vkvk correlators by overwriting the old vkvk.
    To improve a0a0:
    If you use a0p -> a0a0 - 2ca a0p
    If you use pa0 -> a0a0 + 2ca a0p
"""
function improve_corr_vkvk!(vkvk::Vector{uwreal}, vkt0k::Vector{uwreal}, cv::Union{Float64,uwreal})
    der_t0tk = (vkt0k[3:end] - vkt0k[1:end-2]) / 2
    vkvk[2:end-1] = vkvk[2:end-1] .+ 2 .* cv .* der_t0tk
end
improve_corr_vkvk!(vkvk::Corr, t0tk::Corr, cv::Union{Float64,uwreal}) = improve_corr_vkvk!(vkvk.obs, t0tk.obs, cv)

@doc raw"""
ZV(beta::Float64)

Given the coupling beta, this function returns the ZV renormalization constant based on the interpolation formulae
reported in 1811.08209. 
"""
function ZV(beta::Float64)
    g2 = 6/beta
    # CC = [[1.11887e-8, -1.34335e-4*1.11887e-4*1.61288e-4, 4.02275e-5*1.11887e-4*1.44633e-5] [-1.34335e-4*1.11887e-4*1.61288e-4, 1.61288e-8, -4.82986e-5*1.44633e-5*1.61288e-4] [4.02275e-5*1.11887e-4*1.44633e-5, -4.82986e-5*1.44633e-5*1.61288e-4, 1.44633e-10 ]]
    CC = [[1.31619, 4.92750 , 6.15758] [4.92750, 66.8321, 75.3218] [ 6.15758,  75.3218, 85.2733 ]] .*1e-6
    p = cobs([-0.2542, -0.0961, -0.4796], CC, "ZV interp 1")
    zv = 1 - 0.10057*g2 *(1 + p[1]*g2 + p[2]*g2^2)/(1 +p[3]*g2)
    return zv
end

@doc raw"""
cv_loc(beta::Float64)

Given the coupling beta, this function returns the cv improvement coefficient 
for the local vector current based on the results of  1811.08209. 
"""
function cv_loc(beta::Float64)
    
    g2 = 6/beta
    p = uwreal([0.15, 0.35], "cv_loc 1-loop")
    cv_loc = -0.01030 * 4 / 3 * g2 * (1. + p*g2)
    return cv_loc
end
