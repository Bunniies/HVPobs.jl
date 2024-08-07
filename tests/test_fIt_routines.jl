using ADerrors, HVPobs

function myfunc(x, p)
    return p[1] .* exp.(.-x[:,1].*p[2]) .- p[1].*x[:,2]
end

xdata = Array{uwreal}(undef, (20,2))
for k in 1:size(xdata,1)
    xdata[k,1] = uwreal([Float64(k), 0.01*k], "pointx"*string(k))
    xdata[k,2] = uwreal([Float64(k)/2, 0.01*k], "pointxx"*string(k))

end
ydata = Vector{uwreal}(undef, length(xdata[:,1]));

for i in eachindex(ydata)
    ydata[i] = uwreal([myfunc(value.(xdata[i:i,:]), [1.0, 2.0])[1] + 0.01*getindex(randn(1),1), 0.01], "Point "*string(i))
    uwerr(ydata[i]) # We will need the errors for the weights of the fit
end

fit = fit_data(myfunc, xdata, ydata, 2)
uwerr.(fit.param)
fit.param