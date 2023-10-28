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


t00 = comp_t0([ydata,ydata1], [20,80], L=32)


# test rwf deflated
path = "/Users/alessandroconigli/Lattice/data/HVP/rwf_deflated/J303r003.ms1.dat"
read_ms1(path, v="2.0")
get_rw(path,"J303", v="2.0")


#########################
# TEST DATA AUTOMATION
#########################
path = "/Users/alessandroconigli/Lattice/data/HVP/2ptdata"
path_rw = "/Users/alessandroconigli/Lattice/data/aux_obs_data/rwf"
"/Users/alessandroconigli/Lattice/data/aux_obs_data/wilson/H101/H101r000.ms.dat"

cdata = get_data(path, "H101", "light", "V1V1c")
rw = get_rw(path_rw, "H101")
corr = get_corr(path, "H101", "light", "V1V1c")



#########################################
# implement reading_rwf for openQCD 2.0
#########################################
##
path_rw_old = "/Users/alessandroconigli/Lattice/data/aux_obs_data/rwf/N202/N202r001.ms1.dat"
path_rw_new = "/Users/alessandroconigli/Lattice/data/HVP/rwf_deflated/J303r003.ms1.dat"
path_rw_new = "/Users/alessandroconigli/Lattice/data/aux_obs_data/rwf/J500/r006/J500r006.ms1.dat"

tmp_data = HVPobs.Data.read_rw_openQCD2(path_rw_new)
##
data = open(path_rw_new, "r")
nrw = read(data, Int32) # line 120
nrw = Int32(nrw / 2)

nfct = Array{Int32}(undef, nrw)
read!(data, nfct)        # line 138
nfct_inheader = 1

nsrc = Array{Int32}(undef, nrw)
read!(data, nsrc)       # line  146

read(data, Int32)        # 149
# ok up to here  (should be the same as v1.6 apart for nrw = nrw/2)

r_data = Array{Array{Float64}}(undef, nrw)
#[r_data[k] = zeros(Float64, nfct[k], nsrc[k], 1) for k = 1:nrw]

vcfg = Vector{Int32}(undef, 0) 
while !eof(data)
    # read a single config struct
    push!(vcfg, read(data, Int32)) #155
    println("\n ######## cnfg: ", vcfg[end])
    d = read(data, Int32) # 640
    println("d: ", d)
    n = Array{Int32}(undef, d)
    read!(data, n)   #642
    println("n: ", n)
    size = read(data, Int32) # 644
    println("size: ", size)
    
    m = n[1]
    for k in 2:d
        m *= n[k]
    end
    println("m: ", m)

    #tmp_data = Array{Float64}(undef, m)
    #read!(data, tmp_data)

    seek(data,  position(data) + 640 -4 - 4*d -4)
    

end
close(data)

##

glob_head_size = 4 + 4 + 4*nrw*(nfct_inheader + 1)
datsize = 4 + 2*8*sum(nsrc .* nfct)
fsize = filesize(path_rw_new)
# ok up to here

ncnfg = Int32((fsize - glob_head_size)/datsize)



close(data)

# reverse engineer for datsize
datsize_rev = Int32((fsize-glob_head_size)/431)





######################################
# chat GPT Version
########################################

# Define a custom struct for the header data
struct FileHeader
    nrw::Int
    nfct::Vector{Int}
    nsrc::Vector{Int}
end

# Define a custom struct for the data
struct DataStruct
    nc::Int
    sqn::Vector{Vector{Vector{Float64}}}
    lnr::Vector{Vector{Vector{Float64}}}
end

# Function to read the header structure from the file
function read_header(file::IO)
    nrw = read(file, Int)
    nfct = read(file, Vector{Int}, nrw)
    nsrc = read(file, Vector{Int}, nrw)
    return FileHeader(nrw, nfct, nsrc)
end

# Function to read the data structure from the file
function read_data(file::IO)
    nc = read(file, Int)
    sqn = read(file, Vector{Vector{Vector{Float64}}}, file_head.nrw)
    lnr = read(file, Vector{Vector{Vector{Float64}}}, file_head.nrw)
    return DataStruct(nc, sqn, lnr)
end

# Define the file path
file_path = "/Users/alessandroconigli/Lattice/data/HVP/rwf_deflated/J303r003.ms1.dat"  # Update this to your file's path

# Open the file for reading in binary mode
file = open(file_path, "r")

# Read the header data
file_head = read_header(file)

# Read the data for each configuration
data = []
while !eof(file)
    push!(data, read_data(file))
end

# Close the file
close(file)