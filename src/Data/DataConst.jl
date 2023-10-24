const GAMMA = ["V1T10", "V1V1c", "V2T20", "V2V2c", "V3T30", "V3V3c",
               "V1V1", "V1cT10", "V2V2", "V2cT20", "V3V3", "V3cT30",
               "A0P", "PP"  ]


const CLS = Dict(
    "H101" => Dict("repLen" => [1007, 1009], "nms" => 2016),
    "B450" => Dict("repLen" => [1600],       "nms" => 1600),
    "N202" => Dict("repLen" => [899, 1003],  "nms" => 1902),
    "N300" => Dict("repLen" => [1540],       "nms" => 1540),
    "N300" => Dict("repLen" => [1540],       "nms" => 1540),
    "J500" => Dict("repLen" => [789, 655, 431],    "nms" => 1875)
)