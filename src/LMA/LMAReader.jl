function read_eigen_eigen(path::String)
    nn = basename(path)

    f = readdlm(path, skipstart=1)
    # return f

    delim = findall(x-> typeof(x)<:AbstractString && occursin("#", x), f)

    tvals = Int64((size(f)[1] - size(delim)[1]) / size(delim)[1])

    data = Array{Float64}(undef, length(delim), tvals)

    for k in eachindex(delim)
        idx = delim[k].I[1]

        # println(f[idx+1:idx+tvals, 2])
        data[k, :] =  f[idx+1:idx+tvals,2]

    end
    return data
end
