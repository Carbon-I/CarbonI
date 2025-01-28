
using Unitful
using CarbonI

Base.@kwdef mutable struct Scenario
    ### Characteristics of the scenario 
    lat::Float64
    lon::Float64
    sza::Float64
    profile_hr

end


function reference_scenario()
    lat = 0.0
    lon = -62.0
    sza = 30.0
    profile_hr = CarbonI.read_atmos_profile_MERRA2(CarbonI.default_merra_file, lat, lon, 7);
    sce = Scenario(lat=lat, lon=lon, sza=sza, profile_hr=profile_hr)
    return sce
end

function stressing_scenario()
    lat = 0.0
    lon = 0.0
    sza = 0.0
    profile_hr = CarbonI.read_atmos_profile_MERRA2(CarbonI.default_merra_file, lat, lon, 7);
    sce = Scenario(lat=lat, lon=lon, sza=sza, profile_hr=profile_hr)
    return sce
end

function generic_scenario(lat::Float64, lon::Float64, sza::Float64)
    profile_hr = CarbonI.read_atmos_profile_MERRA2(CarbonI.default_merra_file, lat, lon, 7);
    sce = Scenario(lat=lat, lon=lon, sza=sza, profile_hr=profile_hr)
    return sce
end
