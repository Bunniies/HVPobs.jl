using HVPobs
using Revise


path_1cnfg = "/Users/alessandroconigli/Lattice/data/HVP/mesons_charm/B450/Averages/confs_r0/b1/C2pt_n1"
f = HVPobs.Data.read_kappa_charm_data_config(path_1cnfg)

## check kappa values
get_kappa_values()
## check kappa configs
tt = HVPobs.Data.get_kappa_configs()
##
path = "/Users/alessandroconigli/Lattice/data/HVP/mesons_charm/B450"
datatot = read_kappa_charm_all_config(path)

corr_ss = corr_obs(datatot["ss"])
corr_sh1 = corr_obs(datatot["sh1"])
corr_sh4 = corr_obs(datatot["sh4"])

tt = readdir("/Users/alessandroconigli/Lattice/data/HVP/mesons_charm/B450/Averages/confs_r0")