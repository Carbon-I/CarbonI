### A Pluto.jl notebook ###
# v0.20.4

using CarbonI 
include(joinpath(dirname(pathof(CarbonI)), "Requirements", "common.jl"))

using InstrumentOperator, Unitful, Interpolations, DiffResults
using NCDatasets, Polynomials, LinearAlgebra, SpecialPolynomials
using Statistics




scenario = CarbonI.reference_scenario() 
#specs = CarbonI.build_instrument("CBE")
specs = CarbonI.build_instrument("Requirement")

# Create an instrument with the specs
ins = InstrumentOperator.createGratingNoiseModel(specs.ET, specs.Pitch, 
	specs.FPA_quantum_efficiency, specs.bench_efficiency, 
	specs.Fnumber, 2*specs.SSI, 
	(specs.readout_noise), specs.dark_current); 

soil, x, solarIrr, σ_matrix, profile, h, Sₐ = setup_data(scenario, specs)



# Compute baseline spectrally resolved albedo
mean_base_albedo = mean(soil(specs.modelling_wl)[(specs.modelling_wl .>= 2105) .& (specs.modelling_wl .<= 2155)])


# Range of values to sweep
broadband_albedo = 0.04:0.02:0.8
sza_range = 5:5:75

# Create a NetCDF file to store the results
#file_path = "errors_as_function_of_albedo_and_sza_cbe.nc"
file_path = "errors_as_function_of_albedo_and_sza_req.nc"
ds = Dataset(file_path, "c")
# add broadband_albedo to file
dim_albedo = defDim(ds, "albedo", length(broadband_albedo))
dim_sza = defDim(ds, "sza", length(sza_range))

header = ["a","b","c","d","e","f","g","h"]
header_str = join(header, ",")

# Define variables
output_albedo = defVar(ds, "albedo", Float64, ("albedo",))
output_sza = defVar(ds, "sza", Float64, ("sza",))
ch4_error = defVar(ds, "CH4_error", Float64,     ("sza", "albedo"))
co2_error = defVar(ds, "CO2_error", Float64,     ("sza", "albedo"))
h2o_error = defVar(ds, "H2O_error", Float64,     ("sza", "albedo"))
hdo_error = defVar(ds, "HDO_error", Float64,     ("sza", "albedo"))
n2o_error = defVar(ds, "N2O_error", Float64,     ("sza", "albedo"))
co_error = defVar(ds, "CO_error", Float64,       ("sza", "albedo"))
co213_error = defVar(ds, "CO213_error", Float64, ("sza", "albedo"))
c2h6_error = defVar(ds, "C2H6_error", Float64,   ("sza", "albedo"))
output_albedo[:] = broadband_albedo
output_sza[:] = sza_range

output_error = Dict()
for key in ["ch4", "co2", "h2o", "hdo", "n2o", "co", "co213", "c2h6"]
	output_error[key] = zeros(length(sza_range), length(broadband_albedo))
end

# Now loop over scenarioes
alb_ind = 1:length(broadband_albedo)
for (_sza, sza) in enumerate(sza_range)
	Threads.@threads for _alb in alb_ind
		alb = broadband_albedo[_alb]

		# Compute baseline spectrally resolved albedo
		refl   = ones(length(specs.modelling_wl)) * alb;
		error = calc_rel_error(specs, x, solarIrr, refl, sza, σ_matrix, profile, h, ins, Sₐ)

		for key in ["ch4", "co2", "h2o", "hdo", "n2o", "co", "co213", "c2h6"]
			output_error[key][_sza, _alb] = error[key]
		end
		# Compute this for a 400m pixel, as required!
		#rel_ch4_proxy_error_400 = sqrt((n2o_error / sqrt(11.5) / sqrt(400/300) / 330)^2 + (ch4_error / sqrt(11.5) / sqrt(400/300) / 1900)^2) * 100

		@show sza, alb

	end
	#@show sza
end
ch4_error[:] = output_error["ch4"]
co2_error[:] = output_error["co2"]
h2o_error[:] = output_error["h2o"]
hdo_error[:] = output_error["hdo"]
n2o_error[:] = output_error["n2o"]
co_error[:] = output_error["co"]
co213_error[:] = output_error["co213"]
c2h6_error[:] = output_error["c2h6"]
