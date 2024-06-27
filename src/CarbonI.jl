module CarbonI

using vSmartMOM
using vSmartMOM.Absorption
using NCDatasets
using ProgressMeter
#using Plots
using ImageFiltering
using Distributions
using Interpolations
using Polynomials

using vSmartMOM.Absorption
using Parameters
using CSV
using DelimitedFiles
using ForwardDiff, DiffResults

include("tools.jl")
include("forwardModel.jl")
include("GeosChem/GeosChemReader.jl")
include("GeosChem/tools.jl")

end # module CarbonI
