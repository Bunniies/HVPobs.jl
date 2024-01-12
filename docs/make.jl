using Documenter, HVPobs
makedocs(sitename="My Documentation", doctest=true,
        pages=[
               "Home" => "index.md",
               "Data" => "data/index_data.md",
               "Obs"  => "obs/index_obs.md",
               "Fits" => "fits/index_fits.md",
               "API"  => "api/index_api.md"
              ],
        format = Documenter.HTML(prettyurls=false)
       )

push!(LOAD_PATH, "../src/")
