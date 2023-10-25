using Revise
using HVPobs, ADerrors
using PyPlot
using TimerOutputs

path_data = "/Users/alessandroconigli/Lattice/data/HVP/2ptdata"
path_rw = "/Users/alessandroconigli/Lattice/data/aux_obs_data/rwf"

enslist = ["H101"]

@time begin

v1v1 = get_corr(path_data, enslist[1], "light", "V1V1", path_rw)
v2v2 = get_corr(path_data, enslist[1], "light", "V2V2", path_rw)
v3v3 = get_corr(path_data, enslist[1], "light", "V3V3", path_rw)

v1t10 = get_corr(path_data, enslist[1], "light", "V1T10", path_rw)
v2t20 = get_corr(path_data, enslist[1], "light", "V2T20", path_rw)
v3t30 = get_corr(path_data, enslist[1], "light", "V3T30", path_rw)

cv = cv_loc(CLS_db[enslist[1]]["beta"])

improve_corr_vkvk!(v1v1, v1t10, cv)
improve_corr_vkvk!(v2v2, v2t20, cv)
improve_corr_vkvk!(v3v3, v3t30, cv)
end


v1v1_unimp = get_corr(path_data, enslist[1], "light", "V1V1", path_rw)

