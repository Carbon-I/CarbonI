module CarbonI

###### Julia packages to import  ############
using vSmartMOM
using vSmartMOM.Absorption
using NCDatasets
using ProgressMeter
#using Plots
using ImageFiltering
using Distributions
using Interpolations
using Polynomials
using Parameters
using DocStringExtensions      # Documentation
using CSV
using DelimitedFiles
using ForwardDiff, DiffResults
using Artifacts


include("tools.jl")

artifact_file=joinpath(dirname(pathof(CarbonI)), "..", "Artifacts.toml")
xs_folder = artifact_path(artifact_hash("cross_sections", artifact_file))
merra_folder = artifact_path(artifact_hash("merra", artifact_file))
solar_file = joinpath(dirname(pathof(CarbonI)), "..", "data", "solar_irr.nc")
albedo_file = joinpath(dirname(pathof(CarbonI)), "..", "data", "albedo.csv")


include("priorCovariances.jl")
include("instrument_shapes.jl")
include("forwardModel.jl")
include("GeosChem/GeosChemReader.jl")
include("GeosChem/tools.jl")

end # module CarbonI
