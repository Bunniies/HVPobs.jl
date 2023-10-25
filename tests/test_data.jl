using HVPobs
using DelimitedFiles
using ADerrors, PyPlot

##########################
## TESTS WITH 1 replica
##########################

path = "/Users/alessandroconigli/Lattice/data/HVP/2ptdata/N202/light/raw_data/data_HVP_PP"
cdata = read_hvp_data(path, "N202")
# corr no rwf
corr = corr_obs(cdata, rw=nothing)
plot(value.(corr.obs))
display(gcf())
close()
m_av = meff(corr.obs, [20,50], pl=true)

# corr with rwf
path_rwf = "/Users/alessandroconigli/Lattice/data/aux_obs_data/rwf/N202/N202r001.ms1.dat"
rw = read_ms1(path_rwf)
corr_rw = corr_obs(cdata, rw=rw)

@. model(x,p) = p[1] + p[2]*x # p[2]*exp(-p[3]*(x-18))
fit = fit_routine(model, collect(20:40), m_eff[20:40], 2, pval=true)

m_av, m_eff = meff(corr_rw.obs, [20,50], data=true, pl=true)

##########################
## TESTS WITH 2 replicas
##########################

path = "/Users/alessandroconigli/Lattice/data/HVP/2ptdata/H101/light/raw_data/data_HVP_V1V1"

# test reading routine and corr_obs without rwf
cdata = read_hvp_data(path, "H101")
corr = corr_obs(cdata, real=true, rw=nothing)

m_av = meff(corr.obs, [15, 35], pl=true)

# test rwf reading routines and corr_obs with rwf
path_rwf = "/Users/alessandroconigli/Lattice/data/aux_obs_data/rwf/H101/H101r000.ms1.dat"
path_rwf1 = "/Users/alessandroconigli/Lattice/data/aux_obs_data/rwf/H101/H101r001.ms1.dat"
rw = read_ms1(path_rwf)
rw1 = read_ms1(path_rwf1)
corr_rw = corr_obs(cdata, real=true, rw=[rw, rw1])
# test read_ms routine and comp_t0

path_ms = "/Users/alessandroconigli/Lattice/data/aux_obs_data/wilson/H101/H101r000.ms.dat"
path_ms1 = "/Users/alessandroconigli/Lattice/data/aux_obs_data/wilson/H101/H101r001.ms.dat"

ydata = read_ms(path_ms)
ydata1 = read_ms(path_ms1)


t0 = comp_t0([ydata,ydata1], [20,50], L=32)


# test rwf deflated
path = "/Users/alessandroconigli/Lattice/data/HVP/rwf_deflated/J303r003.ms1.dat"
read_ms1(path)


#########################
# TEST DATA AUTOMATION
#########################
path = "/Users/alessandroconigli/Lattice/data/HVP/2ptdata"
path_rw = "/Users/alessandroconigli/Lattice/data/aux_obs_data/rwf"
"/Users/alessandroconigli/Lattice/data/aux_obs_data/wilson/H101/H101r000.ms.dat"

cdata = get_data(path, "H101", "light", "V1V1c")
rw = get_rw(path_rw, "H101")
corr = get_corr(path, "H101", "light", "V1V1c")