### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ ddbfa6eb-6233-4b48-bf37-18023d54fb9d
begin
	using Pkg
	cd("../..")
	Pkg.activate("./")
	using LazyArtifacts
end

# ╔═╡ b71f5796-ff9d-4a0c-b973-d4de8463bdfe
using Revise

# ╔═╡ 19d32bc0-8a70-4f66-a0a5-7054799df57c
# Import CarbonI Source Project here
begin
	using CarbonI
	include(joinpath(dirname(pathof(CarbonI)), "Requirements", "common.jl"))

	using InstrumentOperator, Unitful, Interpolations, DiffResults
	using NCDatasets, Polynomials, LinearAlgebra, SpecialPolynomials
	using Statistics
	using PrettyTables
end 

# ╔═╡ d687e5b4-7bb4-42e0-b150-374e43790254
begin
	using Plots
	using PlutoUI
	plotly(); 
end
	

# ╔═╡ 82179d5c-2c9f-4671-b6f9-e0a61a113527
begin

	# Define Instrument Specifications for Requirement and CBE
	
	req_specs = CarbonI.build_instrument("Requirement") 
	cbe_specs = CarbonI.build_instrument("CBE") 
	
	ins_req = InstrumentOperator.createGratingNoiseModel(req_specs.ET, req_specs.Pitch, 
		req_specs.FPA_quantum_efficiency, req_specs.bench_efficiency, 
		req_specs.Fnumber, 2*req_specs.SSI, 
		(req_specs.readout_noise), req_specs.dark_current); 
	ins_cbe = InstrumentOperator.createGratingNoiseModel(cbe_specs.ET, cbe_specs.Pitch, 
		cbe_specs.FPA_quantum_efficiency, cbe_specs.bench_efficiency, 
		cbe_specs.Fnumber, 2*cbe_specs.SSI, 
		(cbe_specs.readout_noise), cbe_specs.dark_current); 
	



end

# ╔═╡ 7ad80f67-f716-438b-85f3-958aae70120f
begin
	# Set up scenario
	
	#scenario = CarbonI.stressing_scenario() 
	scenario = CarbonI.reference_scenario() 

	soil_req, x_req, solarIrr_req, σ_matrix_req, profile_req, h_req, Sₐ_req = setup_data(scenario, req_specs)
	soil_cbe, x_cbe, solarIrr_cbe, σ_matrix_cbe, profile_cbe, h_cbe, Sₐ_cbe = setup_data(scenario, cbe_specs);

	mean_base_albedo = mean(soil_cbe(cbe_specs.modelling_wl)[(cbe_specs.modelling_wl .>= 2105) .& (cbe_specs.modelling_wl .<= 2155)]);
end

# ╔═╡ 57adfa63-6bb5-42f6-8f8f-f2bc4486fc92
begin
	refl_req   = soil_req(req_specs.modelling_wl) * scenario.broadband_albedo/mean_base_albedo;
	refl_cbe   = soil_cbe(cbe_specs.modelling_wl) * scenario.broadband_albedo/mean_base_albedo
	error_req, F_req = calc_rel_error(req_specs, x_req, solarIrr_req, refl_req, scenario.sza, σ_matrix_req, profile_req, h_req, ins_req, Sₐ_req, return_F=true) 
	error_cbe, F_cbe = calc_rel_error(cbe_specs, x_cbe, solarIrr_cbe, refl_cbe, scenario.sza, σ_matrix_cbe, profile_cbe, h_cbe, ins_cbe, Sₐ_cbe, return_F=true) 
end

# ╔═╡ dde12853-6ea9-4e5e-a413-3e04dce9355a
begin
	# In ppb
	header = ["Species","Albedo", "Wind", "Global Mode CBE", "Global Mode Req", "Target Mode CBE", "Target Mode Req"]
	data_ppb = []
	for key in ["ch4", "co2", "h2o", "hdo", "n2o", "co", "c2h6"]
		push!(data_ppb, [key, scenario.broadband_albedo, 
								scenario.wind_speed, 
								error_cbe[key] / sqrt(cbe_specs.coadd_rate), 
								error_req[key] / sqrt(req_specs.coadd_rate), 
								error_cbe[key], 
								error_req[key]])
	end 
end

# ╔═╡ 8da5074c-5822-4e76-8a1d-d698482bccb0
begin
	# In kg / hr
	data = []
	for key in ["ch4", "co2", "co", "h2o", "hdo", "n2o", "c2h6"]
		push!(data, [key, scenario.broadband_albedo, 
								scenario.wind_speed, 
								CarbonI.jacobs_eq(key, scenario.wind_speed, error_cbe[key] / sqrt(cbe_specs.coadd_rate), cbe_specs.pixel_size_global, style="detect"),
								CarbonI.jacobs_eq(key, scenario.wind_speed, error_req[key] / sqrt(cbe_specs.coadd_rate), cbe_specs.pixel_size_global, style="detect"),
								CarbonI.jacobs_eq(key, scenario.wind_speed, error_cbe[key], req_specs.pixel_size_target, style="detect"),
								CarbonI.jacobs_eq(key, scenario.wind_speed, error_req[key], req_specs.pixel_size_target, style="detect")])
	end 
end

# ╔═╡ a07570ed-c906-4fc3-b66c-493bda8dec5f
data[1]

# ╔═╡ ad838ad4-89a7-4554-856e-135d7b8a3e57
function uts(q; digits=2)
	if typeof(q) == String
		return q
	elseif typeof(q) == Float64
		return round(q, digits=digits)
	else
    	return round(typeof(q), q,digits=digits)
	end
end

# ╔═╡ 11845a8e-2bcd-44f2-9616-d8ec5fa9ab4f
begin
	header_str = " | " * join(header, " | ") * " | ";
	break_str = join([":--- " for x in 1:length(header)], " | ");
	data_str = [[uts(x) for x in d] for d in data];
	data_ppb_str = [[uts(x) for x in d] for d in data_ppb];
end

# ╔═╡ 4a621a49-6a9c-4028-88f3-df5e004d09f4
md"""
| $(header[1]) | $(header[2]) | $(header[3]) | $(header[4]) | $(header[5]) | $(header[6]) | $(header[7]) |
:--- | :--- | :--- | :--- | :--- | :--- | :--- 
$(data_str[1][1]) | $(data_str[1][2]) | $(data_str[1][3]) | $(data_str[1][4]) | $(data_str[1][5]) | $(data_str[1][6]) | $(data_str[1][7])
$(data_str[2][1]) | $(data_str[2][2]) | $(data_str[2][3]) | $(data_str[2][4]) | $(data_str[2][5]) | $(data_str[2][6]) | $(data_str[2][7])
$(data_str[3][1]) | $(data_str[3][2]) | $(data_str[3][3]) | $(data_str[3][4]) | $(data_str[3][5]) | $(data_str[3][6]) | $(data_str[3][7])
"""

# ╔═╡ 9b9e9691-1e63-4f2a-a395-0299cffb8e33
md"""
| $(header[1]) | $(header[2]) | $(header[3]) | $(header[4]) (ppb) | $(header[5]) (ppb) | $(header[6]) (ppb) | $(header[7]) (ppb) |
:--- | :--- | :--- | :--- | :--- | :--- | :--- 
$(data_ppb_str[1][1]) | $(data_ppb_str[1][2]) | $(data_ppb_str[1][3]) | $(data_ppb_str[1][4]) | $(data_ppb_str[1][5]) | $(data_ppb_str[1][6]) | $(data_ppb_str[1][7])
$(data_ppb_str[2][1]) | $(data_ppb_str[2][2]) | $(data_ppb_str[2][3]) | $(data_ppb_str[2][4]) | $(data_ppb_str[2][5]) | $(data_ppb_str[2][6]) | $(data_ppb_str[2][7])
$(data_ppb_str[3][1]) | $(data_ppb_str[3][2]) | $(data_ppb_str[3][3]) | $(data_ppb_str[3][4]) | $(data_ppb_str[3][5]) | $(data_ppb_str[3][6]) | $(data_ppb_str[3][7])
"""

# ╔═╡ db0a1902-6368-452f-ad38-b914fbe14a75
scenario

# ╔═╡ 476345be-cc6c-4725-9f2a-aab9af9dadac
begin
	# Compute noise equivalent radiance:
	nesr_cbe = InstrumentOperator.noise_equivalent_radiance(ins_cbe,(cbe_specs.instrument_wl)u"nm", (F_cbe)u"mW/m^2/nm/sr");
	nesr_cbe_ = nesr_cbe./1u"mW/m^2/nm/sr"

	nesr_req = InstrumentOperator.noise_equivalent_radiance(ins_req,(req_specs.instrument_wl)u"nm", (F_req)u"mW/m^2/nm/sr");
	nesr_req_ = nesr_req./1u"mW/m^2/nm/sr"

	plot(req_specs.instrument_wl, nesr_req_, label="Requirement"); ylabel!("NESR radiance (mW/m²/nm/sr)");
	plot!(cbe_specs.instrument_wl, nesr_cbe_, label="CBE"); ylabel!("NESR radiance (mW/m²/nm/sr)")
	title!("NESR")
end

# ╔═╡ Cell order:
# ╟─ddbfa6eb-6233-4b48-bf37-18023d54fb9d
# ╠═19d32bc0-8a70-4f66-a0a5-7054799df57c
# ╠═d687e5b4-7bb4-42e0-b150-374e43790254
# ╠═b71f5796-ff9d-4a0c-b973-d4de8463bdfe
# ╠═82179d5c-2c9f-4671-b6f9-e0a61a113527
# ╠═7ad80f67-f716-438b-85f3-958aae70120f
# ╠═57adfa63-6bb5-42f6-8f8f-f2bc4486fc92
# ╠═dde12853-6ea9-4e5e-a413-3e04dce9355a
# ╠═8da5074c-5822-4e76-8a1d-d698482bccb0
# ╠═a07570ed-c906-4fc3-b66c-493bda8dec5f
# ╟─ad838ad4-89a7-4554-856e-135d7b8a3e57
# ╟─11845a8e-2bcd-44f2-9616-d8ec5fa9ab4f
# ╟─4a621a49-6a9c-4028-88f3-df5e004d09f4
# ╟─9b9e9691-1e63-4f2a-a395-0299cffb8e33
# ╠═db0a1902-6368-452f-ad38-b914fbe14a75
# ╟─476345be-cc6c-4725-9f2a-aab9af9dadac
