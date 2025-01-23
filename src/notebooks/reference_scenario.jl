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
	include(joinpath(dirname(pathof(CarbonI)), "forwardModel.jl"))
end

# ╔═╡ d687e5b4-7bb4-42e0-b150-374e43790254
begin
	using ImageFiltering, DiffResults, ForwardDiff, InstrumentOperator, Unitful, Interpolations
	using NCDatasets, Polynomials, LinearAlgebra, SpecialPolynomials, DelimitedFiles
	using Plots
	using Artifacts
	using PlutoUI
	plotly(); 
	# Load spectroscopies:
	co2, ch4, h2o, hdo, n2o, co, co2_iso2, c2h6 = CarbonI.loadXSModels();
end
	

# ╔═╡ b298729a-0452-41d2-902d-9be8c0efdc6b
begin
	DS = Dataset(CarbonI.solar_file)
	wlSol = 1e3*DS["wl"][:]
	solar_irr = 1e3*DS["solar_irr"][:] # convert to mW/m2/nm
	close(DS)

	# optionally load a custom reflectance
	#DS = Dataset("data/reflectance_cube_all_1nm_from450nm.h5")
	#r = DS["reflectance_cube"][:]
	#close(DS)
end

# ╔═╡ 750d6e2c-4b3c-45e8-b65e-f7dfb9de2fa8

scenario = CarbonI.reference_scenario()

# ╔═╡ 0a785adf-d7d4-4c78-8b66-217bfdea4ee4
cbe_specs = CarbonI.build_instrument("CBE") 

# ╔═╡ 98395d3e-cfd0-4b3c-b60c-4620f18d7ee3
req_specs = CarbonI.build_instrument("Requirement")

# ╔═╡ 4106bbff-9a3c-4f93-938d-0f32196a2308
begin
	plot(req_specs.modelling_wl .- req_specs.instrument_wl[Int32(round(end/2))], req_specs.convolution_matrix[Int32(round(end/2)),:], label="Requirement")
	
	plot!(cbe_specs.modelling_wl .- cbe_specs.instrument_wl[Int32(round(end/2))], cbe_specs.convolution_matrix[Int32(round(end/2)),:], label="CBE")
	xlims!(-4,4); xlabel!("Δwl (nm)")  
	title!("Instrument Line Shape")
end

# ╔═╡ 4991e247-1e13-4b93-8257-de777824cc95
begin
	hitran_array = (co2, h2o, ch4, co, n2o, hdo, co2_iso2, c2h6);
	# Precompute the cross sections:
	σ_matrix_hr = CarbonI.compute_profile_crossSections(scenario.profile_hr, hitran_array , req_specs.modelling_wl);
	
	nL = length(scenario.profile_hr.T)
		
	vmr_co2 = zeros(nL) .+ 407e-6
	vmr_ch4 = zeros(nL) .+ 1.8e-6
	vmr_ch4[1:3] .= 1.4e-6
	vmr_h2o = scenario.profile_hr.vcd_h2o ./ scenario.profile_hr.vcd_dry
	vmr_co  = zeros(nL) .+ 100e-9
	vmr_n2o = zeros(nL) .+ 337e-9
	vmr_n2o[1:3] .= 100e-9
	vmr_hdo = vmr_h2o * 0.9
	vmr_c2h6 = zeros(nL) .+ 1.0e-9
	vmrs = [vmr_co2, vmr_h2o, vmr_ch4,vmr_co, vmr_n2o, vmr_hdo, vmr_co2, vmr_c2h6];
	
	sol  = CubicSplineInterpolation(range(wlSol[1],wlSol[end], length=length(wlSol)),solar_irr, extrapolation_bc=Interpolations.Flat());
end

# ╔═╡ a43002aa-e378-4b87-ad06-8dacebccabc4
# Reduce dimensions, group layers together to get  layers of equal pressure difference:
n_layers = 10;

# ╔═╡ 7d44f87d-7851-4224-9712-a37a990f4a0d
# Reduce profile and cross sections to fewer dimensions:
profile, σ_matrix, indis, gasProfiles = CarbonI.reduce_profile(n_layers, scenario.profile_hr, σ_matrix_hr,vmrs)

# ╔═╡ 673d3f9b-7a65-4d18-b8e4-b41957e40564
# Number of Legendre Polynomials for the spectrally resolved albedo term
nLeg = 10

# ╔═╡ 68f8bb01-8e44-460f-95fa-4a9f2feb657c
# Define state vector elements for legendre polynomials
xPoly = zeros(nLeg).+eps()

# ╔═╡ 2e9cf224-d29b-432d-ae3e-b2ebd792d957
# First element has to be unity here
xPoly[1] = 1.0

# ╔═╡ baab38dc-8b6a-4f18-9112-4e44c07a9039
# Define full state vector (includes trace gas VMRs and Legendre Polynomial)
x = [reduce(vcat,gasProfiles) ; xPoly ];

# ╔═╡ 173e4607-37de-4b8e-8d2a-63f8f3d35b2f
# Define Covariance matrices and column kernel operator
begin
	# Get prior covariance matrix:
	n_state = length(x);
	Sₐ = zeros(n_state,n_state);
	rel_error = 0.0001;
	# vcd_ratio = profile_caltech.vcd_dry ./ mean(profile_caltech.vcd_dry)
		
	# Fill the diagonal for the trace gases:
	for i=1:80
		Sₐ[i,i] = (rel_error*x[i])^2   
	end
	# CO2 at surface, 100% error
	Sₐ[10,10] = (200x[10])^2
	Sₐ[20,20] = (200x[20])^2
	Sₐ[30,30] = (200x[30])^2
	Sₐ[40,40] = (200x[40])^2
	Sₐ[50,50] = (200x[50])^2
	Sₐ[60,60] = (200x[60])^2
	Sₐ[70,70] = (200x[70])^2
	Sₐ[80,80] = (20000x[80])^2
	# Put in arbitrarily high numbers for the polynomial term, so these won't be constrained at all! 
	for i=81:n_state
		Sₐ[i,i] = 1e2;
	end
	ratio = profile.vcd_dry/sum(profile.vcd_dry);
	
	h_co2 = zeros(length(x));
	h_co213 = zeros(length(x));
	h_ch4 = zeros(length(x));
	h_h2o = zeros(length(x));
	h_co  = zeros(length(x));
	h_hdo = zeros(length(x));
	h_n2o = zeros(length(x));
	h_c2h6 = zeros(length(x));
	
	h_co2[1:10] .= ratio;
	h_h2o[11:20] .= ratio;
	h_ch4[21:30] .= ratio;
	h_co[31:40] .= ratio;
	h_n2o[41:50] .= ratio;
	h_hdo[51:60] .= ratio;
	h_co213[61:70] .= ratio;
	h_c2h6[71:80] .= ratio;
end

# ╔═╡ cabefd7a-1d62-478e-a108-77914078b938
# Load tropical albedo
clima_alb = readdlm(CarbonI.albedo_file,',', skipstart=1);

# ╔═╡ 82179d5c-2c9f-4671-b6f9-e0a61a113527
# Generate albedo function generator
soil = CubicSplineInterpolation(300:2400,clima_alb[:,2]/1.16, extrapolation_bc=Interpolations.Flat());

# ╔═╡ 33f3664a-9423-4293-8c7f-58f5ebc13a30
begin
	# Create solar spectrum at Forward model resolution
	solarIrr_req = sol(req_specs.modelling_wl);
	
	# Compute baseline spectrally resolved albedo
	refl_req   = soil(req_specs.modelling_wl);

	# Create solar spectrum at Forward model resolution
	solarIrr_cbe = sol(cbe_specs.modelling_wl);
	
	# Compute baseline spectrally resolved albedo
	refl_cbe   = soil(cbe_specs.modelling_wl);	
end

# ╔═╡ 1f68d473-7151-4221-a5b3-1ba65381089a
# Baseline Reference albedo to the stressing tropical case
begin
	plot(cbe_specs.instrument_wl, soil(cbe_specs.instrument_wl), label="CBE")
	plot!(req_specs.instrument_wl, soil(req_specs.instrument_wl), label="Requirement")
	title!("Albedo for the stressing tropical case")
end

# ╔═╡ acf49d45-a885-4a5a-8fed-b3bd9990117b
# Run Forward model and Jacobian generation together:
begin
	result_req = DiffResults.JacobianResult(zeros(length(req_specs.instrument_wl)),x);
	ForwardDiff.jacobian!(result_req, (x) -> 
		forward_model_x_(x, 
						 sun=solarIrr_req, 
						 instrument=req_specs.instrument_kernel,
			             reflectance=refl_req, 
			             sza=scenario.sza, 
			             σ_matrix=σ_matrix, 
					     profile=profile,
			             wl=req_specs.modelling_wl), x);
	K_req = DiffResults.jacobian(result_req);
	F_req = DiffResults.value(result_req);  
end

# ╔═╡ 38511b46-3d3a-44a4-a38f-2ee204aeab47
# Run Forward model and Jacobian generation together:
begin
	result_cbe = DiffResults.JacobianResult(zeros(length(cbe_specs.instrument_wl)),x);
	ForwardDiff.jacobian!(result_cbe, (x) -> 
		forward_model_x_(x,  
						 sun=solarIrr_cbe, 
						 instrument=cbe_specs.instrument_kernel,
			             reflectance=refl_cbe, 
			             sza=scenario.sza, 
			             σ_matrix=σ_matrix, 
					     profile=profile,
			             wl=cbe_specs.modelling_wl), x);
	K_cbe = DiffResults.jacobian(result_cbe);
	F_cbe = DiffResults.value(result_cbe); 
end

# ╔═╡ f77a2bdb-a9ae-4f4b-863e-5044b99171f7
begin
	plot(req_specs.instrument_wl, F_req, label="Requirement"); 
	ylabel!("Reflected radiance (mW/m²/nm/sr)") 
	title!("Tropical Scenario")
	plot!(cbe_specs.instrument_wl, F_cbe, label="CBE"); 
end

# ╔═╡ 3f52e140-8383-4504-a9c2-33696f218fe8
# Create an instrument with CBE
cbe_ins = InstrumentOperator.createGratingNoiseModel(cbe_specs.ET, cbe_specs.Pitch, cbe_specs.FPA_quantum_efficiency, cbe_specs.bench_efficiency, cbe_specs.Fnumber, cbe_specs.SSI, (cbe_specs.readout_noise), cbe_specs.dark_current);  

# ╔═╡ e3a94fbb-19a6-42fc-b67c-5f82a0960f6c
# Create an instrument with required parameters
req_ins = InstrumentOperator.createGratingNoiseModel(req_specs.ET, req_specs.Pitch,req_specs.FPA_quantum_efficiency, req_specs.bench_efficiency, req_specs.Fnumber, req_specs.SSI, (req_specs.readout_noise), req_specs.dark_current);

# ╔═╡ 4d07731c-60ee-4f6e-a7fc-76cf129e1bbf
begin
	# Compute noise equivalent radiance:
	nesr_cbe = InstrumentOperator.noise_equivalent_radiance(cbe_ins,(cbe_specs.instrument_wl)u"nm", (F_cbe)u"mW/m^2/nm/sr");
	nesr_cbe_ = nesr_cbe./1u"mW/m^2/nm/sr"

	nesr_req = InstrumentOperator.noise_equivalent_radiance(req_ins,(req_specs.instrument_wl)u"nm", (F_req)u"mW/m^2/nm/sr");
	nesr_req_ = nesr_req./1u"mW/m^2/nm/sr"
end

# ╔═╡ 3eb31d88-8582-43ea-a49c-774a7985f770
begin
	plot(req_specs.instrument_wl, nesr_req_, label="Requirement"); ylabel!("NESR radiance (mW/m²/nm/sr)");
	plot!(cbe_specs.instrument_wl, nesr_cbe_, label="CBE"); ylabel!("NESR radiance (mW/m²/nm/sr)")
	title!("NESR")
end

# ╔═╡ 0e4f27b0-ae5d-4491-ad93-a1eddf8d6110
#Compute electrons at the FPA
#e_cbe = InstrumentOperator.photons_at_fpa(cbe_ins, (cbe_specs.instrument_wl)u"nm", (F_cbe)u"mW/m^2/nm/sr");

# ╔═╡ 7f4d5bf5-3f6e-4e95-a062-2c975eb16535
# Generate S\_epsilon matrix
Se_CBE = Diagonal(nesr_cbe_.^2);

# ╔═╡ 38c196ad-a935-463f-9553-c370683ea3ec
# Generate S\_epsilon matrix
Se_REQ = Diagonal(nesr_req_.^2);

# ╔═╡ 669ede59-e504-48ea-be92-5a1ce7738955
# Compute the Gain Matrix:
G_cbe = inv(K_cbe'inv(Se_CBE)K_cbe + inv(Sₐ))K_cbe'inv(Se_CBE);

# ╔═╡ 18585d46-7917-49df-963f-f8bed76e1435
# Posterior covariance matrix (Profile retrievals, 10 layers):
Ŝ_cbe = inv(K_cbe'inv(Se_CBE)K_cbe + inv(Sₐ));

# ╔═╡ ec5b9496-2616-4f1d-bfa7-89632fe76fcc
# Posterior covariance matrix (Profile retrievals, 10 layers):
Ŝ_req = inv(K_req'inv(Se_REQ)K_req + inv(Sₐ)); 

# ╔═╡ 8c5e3d56-40d7-4666-8f4f-09c78df92b55
begin
	# Apply column operator to get to total column uncertainties:
	ch4_error    = sqrt(h_ch4' * Ŝ_cbe * h_ch4)*1e9
	co2_error    = sqrt(h_co2' * Ŝ_cbe * h_co2)*1e6
	h2o_error    = sqrt(h_h2o' * Ŝ_cbe * h_h2o)*1e6
	hdo_error    = sqrt(h_hdo' * Ŝ_cbe * h_hdo)*1e6
	n2o_error    = sqrt(h_n2o' * Ŝ_cbe * h_n2o)*1e9
	co_error     = sqrt(h_co'  * Ŝ_cbe * h_co)*1e9
	co213_error  = sqrt(h_co213'  * Ŝ_cbe * h_co213)*1e6
	c2h6_error   = sqrt(h_c2h6' * Ŝ_cbe * h_c2h6)*1e9
end

# ╔═╡ 7279e87f-2077-4729-8dae-8e4e52b52412
begin
	# Apply column operator to get to total column uncertainties:
	req_ch4_error    = sqrt(h_ch4' * Ŝ_req * h_ch4)*1e9
	req_co2_error    = sqrt(h_co2' * Ŝ_req * h_co2)*1e6
	req_h2o_error    = sqrt(h_h2o' * Ŝ_req * h_h2o)*1e6
	req_hdo_error    = sqrt(h_hdo' * Ŝ_req * h_hdo)*1e6
	req_n2o_error    = sqrt(h_n2o' * Ŝ_req * h_n2o)*1e9
	req_co_error     = sqrt(h_co'  * Ŝ_req * h_co)*1e9
	req_co213_error  = sqrt(h_co213'  * Ŝ_req * h_co213)*1e6
	req_c2h6_error   = sqrt(h_c2h6' * Ŝ_req * h_c2h6)*1e9
end

# ╔═╡ c6d2d048-3de1-411d-85fe-3fdd902819fa
begin
	# For co-adding:
	@show req_ch4_error/sqrt(10)
	@show req_co2_error/sqrt(10)
	@show req_n2o_error/sqrt(10)
	@show req_n2o_error/sqrt(11.4)/sqrt(400/300) / 330 * 100
	@show n2o_error/sqrt(11.4)
end

# ╔═╡ 411cb18c-3ac3-4a49-9fb3-d55fc9632496
# Compute this for a 400m pixel, as required!

# ╔═╡ 9ee2076c-0de2-4365-b155-b1d58bf0ddc4
rel_ch4_proxy_error_400 = sqrt((req_n2o_error / sqrt(11.5) / sqrt(400/300) / 330)^2 + (req_ch4_error / sqrt(11.5) / sqrt(400/300) / 1900)^2) * 100

# ╔═╡ 65bdc45f-9a7e-4880-af5d-905983f33573
# Expected relatice error at CBE:
rel_ch4_proxy_error_400_CBE = sqrt((n2o_error / sqrt(11.5) / sqrt(400/300) / 330)^2 + (ch4_error / sqrt(11.5) / sqrt(400/300) / 1900)^2) * 100

# ╔═╡ Cell order:
# ╟─ddbfa6eb-6233-4b48-bf37-18023d54fb9d
# ╠═b71f5796-ff9d-4a0c-b973-d4de8463bdfe
# ╠═19d32bc0-8a70-4f66-a0a5-7054799df57c
# ╠═d687e5b4-7bb4-42e0-b150-374e43790254
# ╠═b298729a-0452-41d2-902d-9be8c0efdc6b
# ╠═750d6e2c-4b3c-45e8-b65e-f7dfb9de2fa8
# ╠═0a785adf-d7d4-4c78-8b66-217bfdea4ee4
# ╠═98395d3e-cfd0-4b3c-b60c-4620f18d7ee3
# ╠═4106bbff-9a3c-4f93-938d-0f32196a2308
# ╠═4991e247-1e13-4b93-8257-de777824cc95
# ╠═a43002aa-e378-4b87-ad06-8dacebccabc4
# ╠═7d44f87d-7851-4224-9712-a37a990f4a0d
# ╠═673d3f9b-7a65-4d18-b8e4-b41957e40564
# ╠═68f8bb01-8e44-460f-95fa-4a9f2feb657c
# ╠═2e9cf224-d29b-432d-ae3e-b2ebd792d957
# ╠═baab38dc-8b6a-4f18-9112-4e44c07a9039
# ╟─173e4607-37de-4b8e-8d2a-63f8f3d35b2f
# ╠═cabefd7a-1d62-478e-a108-77914078b938
# ╠═82179d5c-2c9f-4671-b6f9-e0a61a113527
# ╠═33f3664a-9423-4293-8c7f-58f5ebc13a30
# ╠═1f68d473-7151-4221-a5b3-1ba65381089a
# ╠═acf49d45-a885-4a5a-8fed-b3bd9990117b
# ╠═38511b46-3d3a-44a4-a38f-2ee204aeab47
# ╠═f77a2bdb-a9ae-4f4b-863e-5044b99171f7
# ╠═3f52e140-8383-4504-a9c2-33696f218fe8
# ╠═e3a94fbb-19a6-42fc-b67c-5f82a0960f6c
# ╠═4d07731c-60ee-4f6e-a7fc-76cf129e1bbf
# ╠═3eb31d88-8582-43ea-a49c-774a7985f770
# ╠═0e4f27b0-ae5d-4491-ad93-a1eddf8d6110
# ╠═7f4d5bf5-3f6e-4e95-a062-2c975eb16535
# ╠═38c196ad-a935-463f-9553-c370683ea3ec
# ╠═669ede59-e504-48ea-be92-5a1ce7738955
# ╠═18585d46-7917-49df-963f-f8bed76e1435
# ╠═ec5b9496-2616-4f1d-bfa7-89632fe76fcc
# ╠═8c5e3d56-40d7-4666-8f4f-09c78df92b55
# ╠═7279e87f-2077-4729-8dae-8e4e52b52412
# ╠═c6d2d048-3de1-411d-85fe-3fdd902819fa
# ╠═411cb18c-3ac3-4a49-9fb3-d55fc9632496
# ╠═9ee2076c-0de2-4365-b155-b1d58bf0ddc4
# ╠═65bdc45f-9a7e-4880-af5d-905983f33573
