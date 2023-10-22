using HVPobs
using DelimitedFiles

path = "/Users/alessandroconigli/Lattice/data/HVP/H101/light/raw_data/data_HVP_V1V1"



cdata = read_hvp_data(path, "H101")

corr = corr_obs(cdata, real=true, rw=nothing)

m_av = meff(corr.obs, [15, 35], pl=true)

