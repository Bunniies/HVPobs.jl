const GAMMA = ["V1T10", "V1V1c", "V2T20", "V2V2c", "V3T30", "V3V3c",
               "V1V1", "V1cT10", "V2V2", "V2cT20", "V3V3", "V3cT30",
               "A0P", "PP"  ]


const CLS_db = Dict(
    "H101" => Dict("L"=>32, "beta"=>3.4, "kappa_l"=>0.13675962, "kappa_s"=>0.13675962, "dtr"=>2, "plat_t0"=>[20,80]),
    "H102" => Dict("L"=>32, "beta"=>3.4, "kappa_l"=>0.136865, "kappa_s"=>0.136549339, "dtr"=>2, "plat_t0"=>[20,80]),
    "H105" => Dict("L"=>32, "beta"=>3.4, "kappa_l"=>0.136970, "kappa_s"=>0.13634079, "dtr"=>2, "plat_t0"=>[20,80]),
    "N101" => Dict("L"=>48, "beta"=>3.4, "kappa_l"=>0.136970, "kappa_s"=>0.13634079, "dtr"=>2, "plat_t0"=>[20,80]),
    "C101" => Dict("L"=>48, "beta"=>3.4, "kappa_l"=>0.137030, "kappa_s"=>0.136222041, "dtr"=>1, "plat_t0"=>[20,80]),
    
    "B450" => Dict("L"=>32, "beta"=>3.46, "kappa_l"=>0.136890, "kappa_s"=>0.136890, "dtr"=>1, "plat_t0"=>[20,80]),
    "S400" => Dict("L"=>32, "beta"=>3.46, "kappa_l"=>0.136984, "kappa_s"=>0.136702387, "dtr"=>1, "plat_t0"=>[20,80]),
    "N451" => Dict("L"=>48, "beta"=>3.46, "kappa_l"=>0.1370616, "kappa_s"=>0.1365480771, "dtr"=>1, "plat_t0"=>[20,80]),
    "D450" => Dict("L"=>64, "beta"=>3.46, "kappa_l"=>0.137126, "kappa_s"=>0.136420428639937, "dtr"=>1, "plat_t0"=>[25,100]),
    
    "H200" => Dict("L"=>32, "beta"=>3.55, "kappa_l"=>0.137000, "kappa_s"=>0.137000, "dtr"=>1, "plat_t0"=>[20,80]),
    "N202" => Dict("L"=>48, "beta"=>3.55, "kappa_l"=>0.137000, "kappa_s"=>0.137000, "dtr"=>2, "plat_t0"=>[20,110]),
    "N203" => Dict("L"=>48, "beta"=>3.55, "kappa_l"=>0.137080, "kappa_s"=>0.136840284, "dtr"=>1, "plat_t0"=>[20,110]),
    "N200" => Dict("L"=>48, "beta"=>3.55, "kappa_l"=>0.137140, "kappa_s"=>0.13672086, "dtr"=>1, "plat_t0"=>[20,105]),
    "D200" => Dict("L"=>64, "beta"=>3.55, "kappa_l"=>0.137200, "kappa_s"=>0.136601748, "dtr"=>2, "plat_t0"=>[25,100]),
    "E250" => Dict("L"=>96, "beta"=>3.55, "kappa_l"=>0.137232867, "kappa_s"=>0.136536633, "dtr"=>1, "plat_t0"=>[25,170]),
    
    "N300" => Dict("L"=>48, "beta"=>3.7, "kappa_l"=>0.137000, "kappa_s"=>0.137000, "dtr"=>1, "plat_t0"=>[20,105]),
    "N302" => Dict("L"=>48, "beta"=>3.7, "kappa_l"=>0.137064, "kappa_s"=>0.1368721791358, "dtr"=>1, "plat_t0"=>[20,110]),
    "J303" => Dict("L"=>64, "beta"=>3.7, "kappa_l"=>0.137123, "kappa_s"=>0.1367546608, "dtr"=>2, "plat_t0"=>[25,170]),
    "E300" => Dict("L"=>96, "beta"=>3.7, "kappa_l"=>0.137163, "kappa_s"=>0.1366751636177327, "dtr"=>1, "plat_t0"=>[25,170]),
)

const CLS = Dict(
    "H101" => Dict("repLen" => [1007, 1009], "nms" => 2016),
    "B450" => Dict("repLen" => [1600],       "nms" => 1600),
    "N202" => Dict("repLen" => [899, 1003],  "nms" => 1902),
    "N300" => Dict("repLen" => [1540],       "nms" => 1540),
    "N300" => Dict("repLen" => [1540],       "nms" => 1540),
    "J500" => Dict("repLen" => [789, 655, 431],    "nms" => 1875)
)

const b_values = [3.40, 3.46, 3.55, 3.70, 3.85]
const hc = 197.3269804 #MeV fm
# Madrid scale setting
const t0sqrt_ph = uwreal([0.1439, 0.006], "sqrt(t0) [fm]") 
#1608.08900
const t0_data = [2.86, 3.659, 5.164, 8.595, 14.040]
const t0_error = [11, 16, 18, 29, 49] .* 1e-3

const Ct0 = zeros(5, 5)
for i = 1:5
    Ct0[i,i] = t0_error[i] ^ 2    
end

const t0_ = cobs(t0_data, Ct0, "t0sym/a2")
t0(beta::Float64) = t0_[b_values .== beta][1]

