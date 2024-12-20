module CarbonI

using Artifacts, LazyArtifacts
using Pkg.Artifacts: load_artifacts_toml, ensure_artifact_installed
#@show artifact"cross_sections"
#artifact_path(_hash)
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
using CSV
using DelimitedFiles
using ForwardDiff, DiffResults

include("tools.jl")
include("priorCovariances.jl")
include("instrument_shapes.jl")
include("forwardModel.jl")
include("GeosChem/GeosChemReader.jl")
include("GeosChem/tools.jl")

end # module CarbonI
