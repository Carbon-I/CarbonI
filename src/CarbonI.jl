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
using Pkg.Artifacts


include("tools.jl")
const xs_folder = get_artifacts_path("cross_sections")
const merra_folder = get_artifacts_path("merra")
include("priorCovariances.jl")
include("instrument_shapes.jl")
include("forwardModel.jl")
include("GeosChem/GeosChemReader.jl")
include("GeosChem/tools.jl")

end # module CarbonI
