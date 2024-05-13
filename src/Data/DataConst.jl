using ADerrors

const GAMMA = ["V1T10", "V1V1c", "V2T20", "V2V2c", "V3T30", "V3V3c",
               "V1V1", "V1cT10", "V2V2", "V2cT20", "V3V3", "V3cT30",
               "A0P", "PP", "PA0"  ]


const CLS_db = Dict(
    "A653" => Dict("L"=>24, "beta"=>3.34, "kappa_l"=>0.1365716, "kappa_s"=>0.1365716, "dtr"=>1, "plat_t0"=>[15,35]),
    "A654" => Dict("L"=>24, "beta"=>3.34, "kappa_l"=>0.136750, "kappa_s"=>0.136216193, "dtr"=>1, "plat_t0"=>[15,35]),

    "H101" => Dict("L"=>32, "beta"=>3.4, "kappa_l"=>0.13675962, "kappa_s"=>0.13675962, "dtr"=>2, "plat_t0"=>[20,80]),
    "H102" => Dict("L"=>32, "beta"=>3.4, "kappa_l"=>0.136865, "kappa_s"=>0.136549339, "dtr"=>2, "plat_t0"=>[20,80]),
    "H105" => Dict("L"=>32, "beta"=>3.4, "kappa_l"=>0.136970, "kappa_s"=>0.13634079, "dtr"=>2, "plat_t0"=>[20,80]),
    "N101" => Dict("L"=>48, "beta"=>3.4, "kappa_l"=>0.136970, "kappa_s"=>0.13634079, "dtr"=>2, "plat_t0"=>[20,110]),
    "C101" => Dict("L"=>48, "beta"=>3.4, "kappa_l"=>0.137030, "kappa_s"=>0.136222041, "dtr"=>1, "plat_t0"=>[20,80]),
    "C102" => Dict("L"=>48, "beta"=>3.4, "kappa_l"=>0.13705084580022, "kappa_s"=>0.13612906255557, "dtr"=>1, "plat_t0"=>[20,80]),
    "D150" => Dict("L"=>64, "beta"=>3.4, "kappa_l"=>0.137088, "kappa_s"=>0.13610755, "dtr"=>1, "plat_t0"=>[25,110]),
    
    "B450" => Dict("L"=>32, "beta"=>3.46, "kappa_l"=>0.136890, "kappa_s"=>0.136890, "dtr"=>1, "plat_t0"=>[20,50]),
    "S400" => Dict("L"=>32, "beta"=>3.46, "kappa_l"=>0.136984, "kappa_s"=>0.136702387, "dtr"=>1, "plat_t0"=>[20,80]),
    "N451" => Dict("L"=>48, "beta"=>3.46, "kappa_l"=>0.1370616, "kappa_s"=>0.1365480771, "dtr"=>1, "plat_t0"=>[20,80]),
    "D450" => Dict("L"=>64, "beta"=>3.46, "kappa_l"=>0.137126, "kappa_s"=>0.136420428639937, "dtr"=>1, "plat_t0"=>[25,100]),
    "D451" => Dict("L"=>64, "beta"=>3.46, "kappa_l"=>0.137140, "kappa_s"=>0.136337761, "dtr"=>2, "plat_t0"=>[25,100]),
    "D452" => Dict("L"=>64, "beta"=>3.46, "kappa_l"=>0.137163675, "kappa_s"=>0.136345904546, "dtr"=>1, "plat_t0"=>[25,100]),
    
    "H200" => Dict("L"=>32, "beta"=>3.55, "kappa_l"=>0.137000, "kappa_s"=>0.137000, "dtr"=>1, "plat_t0"=>[20,80]),
    "N202" => Dict("L"=>48, "beta"=>3.55, "kappa_l"=>0.137000, "kappa_s"=>0.137000, "dtr"=>2, "plat_t0"=>[20,110]),
    "N203" => Dict("L"=>48, "beta"=>3.55, "kappa_l"=>0.137080, "kappa_s"=>0.136840284, "dtr"=>1, "plat_t0"=>[20,110]),
    "N200" => Dict("L"=>48, "beta"=>3.55, "kappa_l"=>0.137140, "kappa_s"=>0.13672086, "dtr"=>1, "plat_t0"=>[20,105]),
    "D251" => Dict("L"=>64, "beta"=>3.55, "kappa_l"=>0.137140, "kappa_s"=>0.13672086, "dtr"=>1, "plat_t0"=>[25,100]),
    "D200" => Dict("L"=>64, "beta"=>3.55, "kappa_l"=>0.137200, "kappa_s"=>0.136601748, "dtr"=>2, "plat_t0"=>[25,100]),
    "D201" => Dict("L"=>64, "beta"=>3.55, "kappa_l"=>0.1372067, "kappa_s"=>0.136546844, "dtr"=>1, "plat_t0"=>[25,100]),
    "E250" => Dict("L"=>96, "beta"=>3.55, "kappa_l"=>0.137232867, "kappa_s"=>0.136536633, "dtr"=>1, "plat_t0"=>[25,170]),
    
    "N300" => Dict("L"=>48, "beta"=>3.7, "kappa_l"=>0.137000, "kappa_s"=>0.137000, "dtr"=>1, "plat_t0"=>[20,105]),
    "N302" => Dict("L"=>48, "beta"=>3.7, "kappa_l"=>0.137064, "kappa_s"=>0.1368721791358, "dtr"=>1, "plat_t0"=>[20,110]),
    "J303" => Dict("L"=>64, "beta"=>3.7, "kappa_l"=>0.137123, "kappa_s"=>0.1367546608, "dtr"=>2, "plat_t0"=>[25,170]),
    "J304" => Dict("L"=>64, "beta"=>3.7, "kappa_l"=>0.137130, "kappa_s"=>0.1366569203, "dtr"=>1, "plat_t0"=>[25,170]),
    "E300" => Dict("L"=>96, "beta"=>3.7, "kappa_l"=>0.137163, "kappa_s"=>0.1366751636177327, "dtr"=>1, "plat_t0"=>[25,170]),

    "J500" => Dict("L"=>64, "beta"=>3.85, "kappa_l"=>0.136852, "kappa_s"=>0.136852, "dtr"=>2, "plat_t0"=>[25,170]),
    "J501" => Dict("L"=>64, "beta"=>3.85, "kappa_l"=>0.1369032, "kappa_s"=>0.136749715, "dtr"=>1, "plat_t0"=>[25,170])

)
const CLS_kappa_crit = Dict( # taken from 2211.03744
    3.34 => 0.1366938,
    3.4  => 0.1369153,
    3.46 => 0.1370613,
    3.55 => 0.1371715,
    3.7  => 0.1371530,
    3.85 => 0.1369767
)

const CLS_CNFG = Dict(

    "A653" => Dict("repLen" => OrderedDict("r0" => 5050), "nms" => 5050),
    "A654" => Dict("repLen" => OrderedDict("r0" => 5068), "nms" => 5068),

    "H101" => Dict("repLen" => OrderedDict("r0" => 1007, "r1" => 1009), "nms" => 2016),
    "H102" => Dict("repLen" => OrderedDict("r1" => 1029, "r2" => 1008), "nms" => 2037),
    #"H105" => Dict("repLen" => OrderedDict("r1" => 1027, "r2" => 1042), "nms" => 2069),
    "H105" => Dict("repLen" => OrderedDict("r1" => 1027, "r2" => 1042, "r5" => 837), "nms" => 2906),
    "N101" => Dict("repLen" => OrderedDict("r1" => 280, "r3" => 404, "r4" => 240, "r5" => 352, "r6" => 320), "nms" => 1596),
    "C101" => Dict("repLen" => OrderedDict("r14" => 2000, "r15" => 601), "nms" => 2601),
    "C102" => Dict("repLen" => OrderedDict("r3" => 500, "r4" => 500, "r5" => 500), "nms" => 1500),
    "D150" => Dict("repLen" => OrderedDict("r0" => 404), "nms" => 404),

    "B450" => Dict("repLen" => OrderedDict("r0" => 1612), "nms" => 1612),
    "S400" => Dict("repLen" => OrderedDict("r0" => 872, "r1" => 2001), "nms" => 2873),
    "N451" => Dict("repLen" => OrderedDict("r0" => 1011), "nms" => 1011),
    # "D450" => Dict("repLen" => OrderedDict("r10" => 500, "r11" => 1000), "nms" => 1500),
    "D450" => Dict("repLen" => OrderedDict("r10" => 500), "nms" => 500),
    "D451" => Dict("repLen" => OrderedDict("r0" => 1028), "nms" => 1028),
    "D452" => Dict("repLen" => OrderedDict("r1" => 161, "r2" => 1000), "nms" => 1161),

    "N202" => Dict("repLen" => OrderedDict("r1" => 899), "nms" => 899),
    "N203" => Dict("repLen" => OrderedDict("r0" => 756, "r1" => 787),  "nms" => 1543),
    "N200" => Dict("repLen" => OrderedDict("r0" => 856, "r1" => 856),  "nms" => 1712),
    "D251" => Dict("repLen" => OrderedDict("r0" => 403, "r1" => 1610),  "nms" => 2013),
    "D200" => Dict("repLen" => OrderedDict("r0" => 2001),      "nms" => 2001),
    "D201" => Dict("repLen" => OrderedDict("r1" => 1078),      "nms" => 1078),
    "E250" => Dict("repLen" => OrderedDict("r0" => 151, "r1" => 1009), "nms" => 1160),

    "N300" => Dict("repLen" => OrderedDict("r1" => 507, "r2" => 1540), "nms" => 2047),
    "N302" => Dict("repLen" => OrderedDict("r1" => 2201),      "nms" => 2201),
    "J303" => Dict("repLen" => OrderedDict("r3" => 1073),      "nms" => 1073),
    "J304" => Dict("repLen" => OrderedDict("r0" => 830, "r1" => 804),      "nms" => 1634),
    "E300" => Dict("repLen" => OrderedDict("r1" => 1137),      "nms" => 1137),
    
    "J500" => Dict("repLen" => OrderedDict("r4" => 789, "r5" => 655, "r6" => 431),    "nms" => 1875),
    "J501" => Dict("repLen" => OrderedDict("r1" => 1635, "r2" => 1142, "r3" => 1150), "nms" => 3927)
)


