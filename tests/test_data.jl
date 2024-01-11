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
m_av = meff(corr.obs, [15,40], pl=true)

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


t00 = comp_t0([ydata,ydata1], [20,80], L=32, pl=true)


# test rwf deflated
path = "/Users/alessandroconigli/Lattice/data/HVP/rwf_deflated/J303r003.ms1.dat"
read_ms1(path, v="2.0")
get_rw(path,"J303", v="2.0")


#########################
# TEST DATA AUTOMATION
#########################
path = "/Users/alessandroconigli/Lattice/data/HVP/2ptdata"
path_rw = "/Users/alessandroconigli/Lattice/data/HVP/rwf_deflated"
path_ms = "/Users/alessandroconigli/Lattice/data/aux_obs_data/wilson/"

cdata = get_data(path, "H101", "light", "V1V1c")
rw = get_rw(path_rw, "H101", v="")


corr = get_corr(path, EnsInfo("N202"), "light", "V1V1c", path_rw=path_rw)
t0_ens = get_t0(path_ms, "H101")



#########################################
# implement reading_rwf for openQCD 2.0
#########################################
##
path_rw_old = "/Users/alessandroconigli/Lattice/data/aux_obs_data/rwf/N202/N202r001.ms1.dat"
path_rw_new = "/Users/alessandroconigli/Lattice/data/HVP/rwf_deflated/J303r003.ms1.dat"
path_rw_new = "/Users/alessandroconigli/Lattice/data/aux_obs_data/rwf/J500/r006/J500r006.ms1.dat"

tmp_data = HVPobs.Data.read_rw_openQCD2(path_rw_new)
##

############
# test reading E250 (gaps in measurements) 
#########

path_data = "/Users/alessandroconigli/Lattice/data/HVP/2ptdata"
cdat = get_data(path_data, "H101", "light", "V1V1")

path_rw   = "/Users/alessandroconigli/Lattice/data/HVP/rwf_deflated"
rw = get_rw(path_rw, "H101")

## reproduce corr_obs
data = cdat.re_data
rep_len = cdat.rep_len 
vcfg = [cdat.idm[1+sum(rep_len[1:k-1]):sum(rep_len[1:k])] for k in eachindex(rep_len)]
replica = Int64.(maximum.(vcfg))

nms = sum(replica)

idm = cdat.idm[:]
if length(cdat.rep_len) != 1
    idm_sum = [fill((k-1)*sum(replica[1:k-1]), rep_len[k]) for k in eachindex(replica)]
    idm .+= vcat(idm_sum...)
end

tvals = size(data, 2)
obs = [uwreal(data[:,t], cdat.id, replica, idm, nms ) for t in 1:tvals];

data_r, W = HVPobs.Obs.apply_rw(data, rw, cdat.idm, rep_len )
ow = [uwreal(data_r[:,t], cdat.id, replica, idm, nms ) for t in 1:tvals];
W_obs = uwreal(W, cdat.id, replica, idm, nms )

##






corr = corr_obs(cdat, rw=rw)
uwerr.(corr.obs)
mchist(corr.obs[10], "H101")[2000:end]