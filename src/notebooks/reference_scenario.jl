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

# ╔═╡ 5b2c23a7-e8cb-4259-9c9b-67a6fbb439e3
using Revise

# ╔═╡ 19d32bc0-8a70-4f66-a0a5-7054799df57c
# Import CarbonI Source Project here
using CarbonI

# ╔═╡ d687e5b4-7bb4-42e0-b150-374e43790254
begin
	using ImageFiltering, DiffResults, ForwardDiff, InstrumentOperator, Unitful, Interpolations
	using NCDatasets, Polynomials, LinearAlgebra, SpecialPolynomials, DelimitedFiles
	using Plots
	using Artifacts
	# Load spectroscopies:
	co2, ch4, h2o, hdo, n2o, co, co2_iso2, c2h6 = CarbonI.loadXSModels(CarbonI.xs_folder);
end
	

# ╔═╡ 016ec9e4-a7a7-4d62-a666-7bc1d55d36ac
include("src/forwardModel.jl")

# ╔═╡ ca415390-fbb1-4812-8562-1130cae62fc2
plotly();

# ╔═╡ b298729a-0452-41d2-902d-9be8c0efdc6b
begin
	DS = Dataset("data/solar_irr.nc")
	wlSol = 1e3*DS["wl"][:]
	solar_irr = 1e3*DS["solar_irr"][:] # convert to mW/m2/nm
	close(DS)
	
	DS = Dataset("data/reflectance_cube_all_1nm_from450nm.h5")
	r = DS["reflectance_cube"][:]
	close(DS)
end

# ╔═╡ a6da19ca-6bc8-4ac7-8b10-66a041a9cf22
# Load some profile from MERRA:
MD = CarbonI.merra_folder*"/MERRA2_300.tavg3_3d_asm_Nv.20100610.nc4"

# ╔═╡ cf06e4b3-4149-422a-9dc9-76a76e2d0b37
# Choose spot in the Amazon (latitude)
myLat = 0.0;

# ╔═╡ ff03de75-d086-4eb1-9288-cc4ac30cbd6f
# Choose spot in the Amazon (longitude)
myLon = -62

# ╔═╡ 2188edda-48c0-43d2-90d9-89b887cdad9f
# Read in high resolution profile for the tropics:
profile_hr = CarbonI.read_atmos_profile_MERRA2(MD, myLat, myLon, 7);

# ╔═╡ 5842480f-50ad-484e-b3c6-d8f5c82a2e70
# Define high resolution wavelength spacing (in nm)
Δwl = 0.005;

# ╔═╡ 53b88229-5974-46aa-964f-4abf57673738
# Define high resolution wavelength array for the RT forward model
wl = 2035:Δwl:2385;

# ╔═╡ d3fe98ca-5996-4332-aabb-489ced362fa2
# Instrument Line Shape as Required
# Create convolution matrix cM and CarbonI sampling grid wl_ci_req
cM_req, wl_ci = CarbonI.create_carbonI_conv_matrix(wl);

# ╔═╡ 64d02acd-3061-420f-9b46-b151e0052cbc
# Instrument Line Shape as CBE
# Create convolution matrix cM and CarbonI sampling grid wl_ci_req
cM_cbe, _ = CarbonI.create_carbonI_conv_matrix_cbe(wl);

# ╔═╡ 5535969a-ce85-41b8-bb56-722cd67db666
begin
	# Define an instrument:
	# LSF(λ)
	FWHM  = 2.2  # FWHM  = 2.2 creates an effective FWHM of 2.5nm, i.e. required. (1.5 creates an effective kernel of 1.88)
	# Spectral Sampling Interval
	SSI_  = 0.7
	# Slit blur (2*SSI box)
	kern1 = CarbonI.box_kernel(2*SSI_, Δwl)
	# LSF (Gaussian)
	kern2 = CarbonI.gaussian_kernel(FWHM, Δwl)
	# Pixel response (1*SSI)
	kern3 = CarbonI.box_kernel(SSI_, Δwl)
	
	# Combine the kernels:
	kernf = imfilter(imfilter(kern1, kern2), kern3)
	# Create the instrument response:
	lociBox = CarbonI.KernelInstrument(kernf, collect(2040:SSI_:2380));
end

# ╔═╡ 4106bbff-9a3c-4f93-938d-0f32196a2308
begin
	plot(wl .- wl_ci[250], cM_req[250,:], label="Required ILS")
	plot!(wl .- wl_ci[250], cM_cbe[250,:], label="CBE ILS")
	xlims!(-4,4); xlabel!("Δwl (nm)")
end

# ╔═╡ 1fe9f8ac-8442-4e64-8f62-b808553ff20b
begin
	# Define the fixed (design choices) for the instrument specs at 400m:
	ET  = 44.0u"ms"         # Exposure time
	SSI = (2*0.7)u"nm"      # Spectral resolution
	Pitch = 18.0u"μm"       # Pixel pitch
	Fnumber = 2.2           # F-number
end

# ╔═╡ 3b2e02d0-3199-41e0-bf8a-8750c72bc300
# CBE FPA quantum efficiency
cbe_FPA_QE = 0.85;           

# ╔═╡ f28267b7-e910-4396-b9e0-baa9fc7896db
# required FPA quantum efficiency
req_FPA_QE = 0.80;           

# ╔═╡ 8ace0941-500e-44f3-ae5f-7386e4f6b638
# CBE Bench efficiency
cbe_Bench_efficiency = 0.72; 

# ╔═╡ a5b11d67-e72f-4888-b256-6986bc571833
# required Bench efficiency
req_Bench_efficiency = 0.72; 

# ╔═╡ 26fbf305-701b-44ad-b297-91979f6bb46b
# CBE Readout noise
cbe_readout_noise = 100;    

# ╔═╡ 6c7dbb18-ce72-4cec-b351-37dcb31a97f5
# required Readout noise
req_readout_noise = 120;    

# ╔═╡ cdf91ae6-1197-4328-bbf9-16046d28c05d
# CBE Dark current
cbe_dark_current = 5e3u"1/s"; 

# ╔═╡ f9dfd67f-243f-4fc0-b992-4dec962ec668
# required Dark current
req_dark_current = 10e3u"1/s"; 

# ╔═╡ 4991e247-1e13-4b93-8257-de777824cc95
begin
	hitran_array = (co2, h2o, ch4, co, n2o, hdo, co2_iso2, c2h6);
	# Precompute the cross sections:
	σ_matrix_hr = CarbonI.compute_profile_crossSections(profile_hr, hitran_array , wl);
	
	nL = length(profile_hr.T)
		
	vmr_co2 = zeros(nL) .+ 407e-6
	vmr_ch4 = zeros(nL) .+ 1.8e-6
	vmr_ch4[1:3] .= 1.4e-6
	vmr_h2o = profile_hr.vcd_h2o ./ profile_hr.vcd_dry
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
profile, σ_matrix, indis, gasProfiles = CarbonI.reduce_profile(n_layers, profile_hr, σ_matrix_hr,vmrs);

# ╔═╡ 4580f545-fc13-4b3a-9ceb-4fd50a8d3451
#size(σ_matrix)

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

# ╔═╡ 9117a48c-1a02-4b2a-a8c8-7fb80d981dbf
# Solar Zenith Angle
sza = 30

# ╔═╡ cabefd7a-1d62-478e-a108-77914078b938
# Load tropical albedo
clima_alb = readdlm("data/albedo.csv",',', skipstart=1);

# ╔═╡ 82179d5c-2c9f-4671-b6f9-e0a61a113527
# Generate albedo function generator
soil = CubicSplineInterpolation(300:2400,clima_alb[:,2]/1.16, extrapolation_bc=Interpolations.Flat());

# ╔═╡ 33f3664a-9423-4293-8c7f-58f5ebc13a30
# Create solar spectrum at Forward model resolution
solarIrr = sol(wl);

# ╔═╡ 21f70ce2-6764-4125-96ea-900f1c106844
# Compute baseline spectrally resolved albedo
refl   = soil(wl);

# ╔═╡ 1f68d473-7151-4221-a5b3-1ba65381089a
# Baseline Reference albedo to the stressing tropical case
plot(wl_ci, soil(wl_ci), label="Albedo for the stressing tropical case")

# ╔═╡ f431645a-cbf0-43fd-a70b-2dc848f3e2e0
result = DiffResults.JacobianResult(zeros(length(lociBox.ν_out)),x);

# ╔═╡ acf49d45-a885-4a5a-8fed-b3bd9990117b
# Run Forward model and Jacobian generation together:
ForwardDiff.jacobian!(result, forward_model_x_, x);K = DiffResults.jacobian(result);F = DiffResults.value(result);

# ╔═╡ f77a2bdb-a9ae-4f4b-863e-5044b99171f7
plot(lociBox.ν_out, F, label="Tropical Scenario"); ylabel!("Reflected radiance (mW/m²/nm/sr)")

# ╔═╡ 3f52e140-8383-4504-a9c2-33696f218fe8
# Create an instrument with CBE
cbe_ins = InstrumentOperator.createGratingNoiseModel(ET, Pitch,cbe_FPA_QE, cbe_Bench_efficiency, Fnumber, SSI, (cbe_readout_noise), cbe_dark_current);

# ╔═╡ e3a94fbb-19a6-42fc-b67c-5f82a0960f6c
# Create an instrument with required parameters
req_ins = InstrumentOperator.createGratingNoiseModel(ET, Pitch,req_FPA_QE, req_Bench_efficiency, Fnumber, SSI, (req_readout_noise), req_dark_current);

# ╔═╡ 4d07731c-60ee-4f6e-a7fc-76cf129e1bbf
# Compute noise equivalent radiance:
nesr_cbe = InstrumentOperator.noise_equivalent_radiance(cbe_ins, (lociBox.ν_out)u"nm", (F)u"mW/m^2/nm/sr");nesr_cbe_ = nesr_cbe./1u"mW/m^2/nm/sr"

# ╔═╡ c5b7d02d-8108-4f60-8f0f-a3e38893e22c
# Compute noise equivalent radiance:
nesr_req = InstrumentOperator.noise_equivalent_radiance(req_ins, (lociBox.ν_out)u"nm", (F)u"mW/m^2/nm/sr");nesr_req_ = nesr_req./1u"mW/m^2/nm/sr"

# ╔═╡ 3eb31d88-8582-43ea-a49c-774a7985f770
begin
	plot(lociBox.ν_out, nesr_req_, label="Required NESR"); ylabel!("NESR radiance (mW/m²/nm/sr)");
	plot!(lociBox.ν_out, nesr_cbe_, label="CBE NESR"); ylabel!("NESR radiance (mW/m²/nm/sr)")
end

# ╔═╡ 0e4f27b0-ae5d-4491-ad93-a1eddf8d6110
#Compute electrons at the FPA
e = InstrumentOperator.photons_at_fpa(cbe_ins, (lociBox.ν_out)u"nm", (F)u"mW/m^2/nm/sr");

# ╔═╡ 7f4d5bf5-3f6e-4e95-a062-2c975eb16535
# Generate S\_epsilon matrix
Se_CBE = Diagonal(nesr_cbe_.^2);

# ╔═╡ 38c196ad-a935-463f-9553-c370683ea3ec
# Generate S\_epsilon matrix
Se_REQ = Diagonal(nesr_req_.^2);

# ╔═╡ 669ede59-e504-48ea-be92-5a1ce7738955
# Compute the Gain Matrix:
G = inv(K'inv(Se_CBE)K + inv(Sₐ))K'inv(Se_CBE);

# ╔═╡ 18585d46-7917-49df-963f-f8bed76e1435
# Posterior covariance matrix (Profile retrievals, 10 layers):
Ŝ_cbe = inv(K'inv(Se_CBE)K + inv(Sₐ));

# ╔═╡ ec5b9496-2616-4f1d-bfa7-89632fe76fcc
# Posterior covariance matrix (Profile retrievals, 10 layers):
Ŝ_req = inv(K'inv(Se_REQ)K + inv(Sₐ));

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
# ╠═19d32bc0-8a70-4f66-a0a5-7054799df57c
# ╟─d687e5b4-7bb4-42e0-b150-374e43790254
# ╟─5b2c23a7-e8cb-4259-9c9b-67a6fbb439e3
# ╠═ca415390-fbb1-4812-8562-1130cae62fc2
# ╠═b298729a-0452-41d2-902d-9be8c0efdc6b
# ╠═016ec9e4-a7a7-4d62-a666-7bc1d55d36ac
# ╠═a6da19ca-6bc8-4ac7-8b10-66a041a9cf22
# ╠═cf06e4b3-4149-422a-9dc9-76a76e2d0b37
# ╠═ff03de75-d086-4eb1-9288-cc4ac30cbd6f
# ╠═2188edda-48c0-43d2-90d9-89b887cdad9f
# ╠═5842480f-50ad-484e-b3c6-d8f5c82a2e70
# ╠═53b88229-5974-46aa-964f-4abf57673738
# ╠═d3fe98ca-5996-4332-aabb-489ced362fa2
# ╠═64d02acd-3061-420f-9b46-b151e0052cbc
# ╠═5535969a-ce85-41b8-bb56-722cd67db666
# ╠═4106bbff-9a3c-4f93-938d-0f32196a2308
# ╠═1fe9f8ac-8442-4e64-8f62-b808553ff20b
# ╠═3b2e02d0-3199-41e0-bf8a-8750c72bc300
# ╠═f28267b7-e910-4396-b9e0-baa9fc7896db
# ╠═8ace0941-500e-44f3-ae5f-7386e4f6b638
# ╠═a5b11d67-e72f-4888-b256-6986bc571833
# ╠═26fbf305-701b-44ad-b297-91979f6bb46b
# ╠═6c7dbb18-ce72-4cec-b351-37dcb31a97f5
# ╠═cdf91ae6-1197-4328-bbf9-16046d28c05d
# ╠═f9dfd67f-243f-4fc0-b992-4dec962ec668
# ╟─4991e247-1e13-4b93-8257-de777824cc95
# ╠═a43002aa-e378-4b87-ad06-8dacebccabc4
# ╠═7d44f87d-7851-4224-9712-a37a990f4a0d
# ╠═4580f545-fc13-4b3a-9ceb-4fd50a8d3451
# ╠═673d3f9b-7a65-4d18-b8e4-b41957e40564
# ╠═68f8bb01-8e44-460f-95fa-4a9f2feb657c
# ╠═2e9cf224-d29b-432d-ae3e-b2ebd792d957
# ╠═baab38dc-8b6a-4f18-9112-4e44c07a9039
# ╟─173e4607-37de-4b8e-8d2a-63f8f3d35b2f
# ╠═9117a48c-1a02-4b2a-a8c8-7fb80d981dbf
# ╠═cabefd7a-1d62-478e-a108-77914078b938
# ╠═82179d5c-2c9f-4671-b6f9-e0a61a113527
# ╠═33f3664a-9423-4293-8c7f-58f5ebc13a30
# ╠═21f70ce2-6764-4125-96ea-900f1c106844
# ╠═1f68d473-7151-4221-a5b3-1ba65381089a
# ╠═f431645a-cbf0-43fd-a70b-2dc848f3e2e0
# ╠═acf49d45-a885-4a5a-8fed-b3bd9990117b
# ╠═f77a2bdb-a9ae-4f4b-863e-5044b99171f7
# ╠═3f52e140-8383-4504-a9c2-33696f218fe8
# ╠═e3a94fbb-19a6-42fc-b67c-5f82a0960f6c
# ╠═4d07731c-60ee-4f6e-a7fc-76cf129e1bbf
# ╠═c5b7d02d-8108-4f60-8f0f-a3e38893e22c
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
