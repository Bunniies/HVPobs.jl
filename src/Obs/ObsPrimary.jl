function corr_obs(cd::CData; real::Bool=true, rw::Union{Array{Float64,2}, Vector{Array{Float64,2}}, Nothing}=nothing, L::Int64=1, nms::Union{Int64, Nothing}=nothing)
    
    real ? data = cd.re_data ./ L^3 : data = cd.im_data ./ L^3
    tvals = size(data, 2)
    nms = isnothing(nms) ?  sum(cd.rep_len) : nms
    
    idm = copy(cd.idm)
    if length(cd.rep_len) != 1
        idm_sum = [fill((k-1)*sum(cd.rep_len[1:k-1]), cd.rep_len[k]) for k in eachindex(cd.rep_len)]
        idm .+= vcat(idm_sum...)
    end

    if isnothing(rw)
        obs = [uwreal(data[:,t], cd.id, cd.rep_len, idm, nms) for t in 1:tvals]
    else
        data_r, W = apply_rw(data, rw, cd.idm)
        ow = [uwreal(data_r[:,t], cd.id, cd.rep_len, idm, nms) for t in 1:tvals]
        W_obs = uwreal(W, cd.id, cd.rep_len, idm, nms)
        obs = [ow[t] / W_obs for t in 1:tvals]
    end

    return Corr(obs, cd)
end

@doc raw"""
    comp_t0(Y::YData, plat::Vector{Int64}; L::Int64, pl::Bool=false, rw::Union{Matrix{Float64}, Nothing}=nothing, npol::Int64=2, wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing, info::Bool=false)

    comp_t0(Y::Vector{YData}, plat::Vector{Int64}; L::Int64, pl::Bool=false, rw::Union{Vector{Matrix{Float64}}, Nothing}=nothing, npol::Int64=2, wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing, info::Bool=false)

Computes `t0` using the energy density of the action `Ysl`(Yang-Mills action).
`t0` is computed in the plateau `plat`.
A polynomial interpolation in `t` is performed to find `t0`, where `npol` is the degree of the polynomial (linear fit by default)

The flag `pl` allows to show the plot. 

The flag `info` provides extra output that contains information about the primary observables. The function returns the primary observables ``<WY>`` and ``<W>``
(it returns the observable <Y> if rw=nothing)

```@example
#Single replica
Y = read_ms(path)
rw = read_ms(path_rw)

t0, Yobs = comp_t0(Y, [38, 58], L=32, info=true)
t0_r, WYobs, Wobs = comp_t0(Y, [38, 58], L=32, rw=rw, info=true)

#Two replicas
Y1 = read_ms(path1)
Y2 = read_ms(path2)
rw1 = read_ms(path_rw1)
rw2 = read_ms(path_rw2)

t0 = comp_t0([Y1, Y2], [38, 58], L=32, pl=true)
t0_r = comp_t0(Y, [38, 58], L=32, rw=[rw1, rw2], pl=true)

```
"""
function comp_t0(Y::YData, plat::Vector{Int64}; L::Int64, pl::Bool=false, 
    rw::Union{Matrix{Float64}, Nothing}=nothing, npol::Int64=2, ws::ADerrors.wspace=ADerrors.wsg,
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing, info::Bool=false)

    Ysl = Y.obs
    t = Y.t
    id = Y.id
    replica = size.([Ysl], 1)

    #Truncation
    if id in keys(ADerrors.wsg.str2id)
        n_ws = findfirst(x-> x == ws.str2id[id], ws.map_nob)
        if !isnothing(n_ws)
            ivrep_ws = ws.fluc[n_ws].ivrep

            if length(ivrep_ws) != 1
                error("Different number of replicas")
            end

            if replica[1] > ivrep_ws[1]
                println("Automatic truncation in Ysl ", ivrep_ws[1], " / ", replica[1], ". R = 1")
                Ysl = Ysl[1:ivrep_ws[1], :, :]
            elseif replica[1] < ivrep_ws[1]
                println(replica[1])
                error("Automatic truncation failed. R = 1\nTry using truncate_data!")
            end
        end
    end

    xmax = size(Ysl, 2)
    nt0 = t0_guess(t, Ysl, plat, L)
    
    dt0 = iseven(npol) ? Int64(npol / 2) : Int64((npol+1)/ 2)
    Y_aux = Matrix{uwreal}(undef, xmax, 2*dt0+1)
    
    if !isnothing(rw)
        Ysl_r, W = apply_rw_t0(Ysl, rw)
        W_obs = uwreal(W, id)
        WY_aux = Matrix{uwreal}(undef, xmax, 2*dt0+1)
    end

    for i = 1:xmax
        k = 1
        for j = nt0-dt0:nt0+dt0
            if isnothing(rw)
                Y_aux[i, k] = uwreal(Ysl[:, i, j], id)
            else
                WY_aux[i, k] = uwreal(Ysl_r[:, i, j], id)
                Y_aux[i, k] = WY_aux[i, k] / W_obs
            end
            k = k + 1
        end
    end
    x = t[nt0-dt0:nt0+dt0]
    t2E = [plat_av(Y_aux[:, j], plat, wpm=wpm) for j=1:2*dt0+1] .* x.^2 / L^3
    
    model(x, p) = get_model(x, p, npol)

    fit  = fit_routine(model, x, t2E, npol)
    par = fit.param
    fmin(x, p) = model(x, p) .- 0.3
    t0 = root_error(fmin, t[nt0], par)
    if pl
        if isnothing(wpm)
            uwerr(t0)
            uwerr.(t2E)
        else
            uwerr(t0, wpm)
            [uwerr(t2E_aux, wpm) for t2E_aux in t2E]
        end

        v = value.(t2E)
        e = err.(t2E)

        figure()
        errorbar(x, v, e, fmt="d", capsize=2, mfc="none", color="black")
        errorbar(value(t0), 0.3, xerr=err(t0), capsize=2, color="red", fmt="d")
        xarr = Float64.(range(5.0, 5.2, length=100))
        yarr = model(xarr, par); uwerr.(yarr)
        fill_between(xarr, value.(yarr) .- err.(yarr), value.(yarr) .+ err.(yarr), color="royalblue", alpha=0.2 )
        ylabel(L"$t^2\langle E(x_0, t)\rangle$")
        xlabel(L"$t/a^2$")
        tight_layout()
        display(gcf())

        figure()
        t_pl = dt0 + 1
        yyy = Y_aux[:, t_pl] * t[nt0]^2 / L^3 ; uwerr.(yyy)
        # plot(collect(1:xmax), value.(Y_aux[:, t_pl]) * t[nt0]^2 / L^3, "x") 
        errorbar(collect(1:xmax), value.(yyy), err.(yyy), fmt="d", mfc="none", color="black", capsize=2) 
        fill_between(collect(plat[1]:plat[2]), v[t_pl]+e[t_pl], v[t_pl]-e[t_pl], alpha=0.75, color="royalblue")
        ylabel(L"$t^2\langle E(x_0, t)\rangle$")
        xlabel(L"$x_0/a$")
        ylim(v[t_pl]-30*e[t_pl], v[t_pl]+50*e[t_pl])
        title(string(L"$t/a^2 = $", t[nt0]))
        tight_layout()
        display(gcf())
        close("all")
    end
    if info && !isnothing(rw)
        return (t0, WY_aux, W_obs)
    elseif info && isnothing(rw)
        return (t0, Y_aux)
    else
        return t0
    end
end

function comp_t0(Y::Vector{YData}, plat::Vector{Int64}; L::Int64, pl::Bool=false, 
    rw::Union{Vector{Matrix{Float64}}, Nothing}=nothing, npol::Int64=2, ws::ADerrors.wspace=ADerrors.wsg,
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing, info::Bool=false)

    nr = length(Y)
    Ysl = getfield.(Y, :obs)
    t = getfield.(Y, :t)
    t = t[1]
    id = getfield.(Y, :id)
    replica = size.(Ysl, 1)

    if !all(id .== id[1])
        error("IDs are not equal")
    end
    id = id[1]
    #Truncation
    if id in keys(ADerrors.wsg.str2id)
        n_ws = findfirst(x-> x == ws.str2id[id], ws.map_nob)
        if !isnothing(n_ws)
            ivrep_ws = ws.fluc[n_ws].ivrep

            if length(replica) != length(ivrep_ws)
                error("Different number of replicas")
            end

            for k in eachindex(replica)
                if replica[k] > ivrep_ws[k]
                    println("Automatic truncation in Ysl ", ivrep_ws[k], " / ", replica[k], ". R = ", k)
                    Ysl[k] = Ysl[k][1:ivrep_ws[k], :, :]
                elseif replica[k] < ivrep_ws[k]
                    error("Automatic truncation failed. R = ", replica[k], "\nTry using truncate_data!")
                end
            end
            replica = size.(Ysl, 1)
        end
    end

    tmp = Ysl[1]
    [tmp = cat(tmp, Ysl[k], dims=1) for k=2:nr]
    nt0 = t0_guess(t, tmp, plat, L)
    xmax = size(tmp, 2)

    dt0 = iseven(npol) ? Int64(npol / 2) : Int64((npol+1) / 2)
    Y_aux = Matrix{uwreal}(undef, xmax, 2*dt0+1)

    if !isnothing(rw)
        Ysl_r, W = apply_rw_t0(Ysl, rw)
        tmp_r = Ysl_r[1]
        tmp_W = W[1]
        [tmp_r = cat(tmp_r, Ysl_r[k], dims=1) for k = 2:nr]
        [tmp_W = cat(tmp_W, W[k], dims=1) for k = 2:nr]
        W_obs = uwreal(tmp_W, id, replica)
        WY_aux = Matrix{uwreal}(undef, xmax, 2*dt0+1)
    end
    for i = 1:xmax
        k = 1
        for j = nt0-dt0:nt0+dt0
            if isnothing(rw)
                Y_aux[i, k] = uwreal(tmp[:, i, j], id, replica)
            else
                WY_aux[i, k] = uwreal(tmp_r[:, i, j], id, replica)
                Y_aux[i, k] = WY_aux[i, k] / W_obs
            end
            k = k + 1
        end
    end
    x = t[nt0-dt0:nt0+dt0]
    t2E = [plat_av(Y_aux[:, j], plat, wpm=wpm) for j=1:2*dt0+1] .* x.^2 / L^3
    
    model(x, p) = get_model(x, p, npol)

    fit = fit_routine(model, x, t2E, npol)
    par = fit.param
    fmin(x, p) = model(x, p) .- 0.3
    t0 = root_error(fmin, t[nt0], par)
    if pl
        if isnothing(wpm)
            uwerr(t0)
            uwerr.(t2E)
        else
            uwerr(t0, wpm)
            [uwerr(t2E_aux, wpm) for t2E_aux in t2E]
        end

        v = value.(t2E)
        e = err.(t2E)

        figure()
        errorbar(x, v, e, fmt="x")
        
        errorbar(value(t0), 0.3, xerr=err(t0), fmt="x")
        ylabel(L"$t^2E$")
        xlabel(L"$t/a^2$")
        display(gcf())

        figure()
        t_pl = dt0 + 1
        plot(collect(1:xmax), value.(Y_aux[:, t_pl]) * t[nt0]^2 / L^3, "x")
        fill_between(collect(plat[1]:plat[2]), v[t_pl]+e[t_pl], v[t_pl]-e[t_pl], alpha=0.75, color="green")
        ylabel(L"$t^2E(x0, t)$")
        xlabel(L"$x_0/a$")
        title(string(L"$t/a^2 = $", t[nt0]))
        display(gcf())
        close("all")
    end
    if info && !isnothing(rw)
        return (t0, WY_aux, W_obs)
    elseif info && isnothing(rw)
        return (t0, Y_aux)
    else
        return t0
    end
end

