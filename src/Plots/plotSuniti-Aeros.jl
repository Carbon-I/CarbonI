using CUDA
device!(1)
using DelimitedFiles
using Distributions
using ImageFiltering
using InstrumentOperator
using Interpolations
using JLD2
using LaTeXStrings
using Plots
using Revise
using Statistics
using vSmartMOM
#using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using vSmartMOM.SolarModel
# Benchmarks: http://web.gps.caltech.edu/~vijay/Rayleigh_Scattering_Tables/STOKES/
SSI = 0.7 #"nm"
FWHM = 1.93
FWHM_gaussian = 0.65 #, # Instead of 0.65!
lower_wavelength = 2036.0
upper_wavelength = 2372.0
##########################
##### Solar spectrum #####
##########################
#=Tsolar = solar_transmission_from_file("/home/sanghavi/code/github/vSmartMOM.jl/src/SolarModel/solar.out")
Tsolar_interp = LinearInterpolation(Tsolar[4:end, 1], Tsolar[4:end, 2])
T_sun = 5777. # K
ν = collect(1e7 ./ (lower_wavelength:SSI:upper_wavelength))  # cm⁻¹
P = planck_spectrum_wn(T_sun, ν) * 2.1629e-05 * π  # mW/m²-cm⁻¹
F₀ = zeros(length(P));
F₀ = Tsolar_interp.(ν) .* P;
=#

parameters =
    parameters_from_yaml("/home/sanghavi/code/github/vSmartMOM.jl/test/test_parameters/aerosol_parameters2.yaml");
model      = model_from_parameters(#fwd_mode,
    parameters);
n_bands = length(parameters.spec_bands)
n_aer = isnothing(parameters.scattering_params) ? 0 : length(parameters.scattering_params.rt_aerosols)
truncation_type = vSmartMOM.Scattering.δBGE{parameters.float_type}(parameters.l_trunc, parameters.Δ_angle)
vmr = isnothing(parameters.absorption_params) ? Dict() : parameters.absorption_params.vmr
p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr =
    vSmartMOM.CoreRT.compute_atmos_profile_fields(parameters.T, parameters.p, parameters.q, vmr)
profile = vSmartMOM.CoreRT.AtmosphericProfile(parameters.T, p_full, parameters.q, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr)
# Reduce the profile to the number of target layers (if specified)
if parameters.profile_reduction_n != -1
    profile = vSmartMOM.CoreRT.reduce_profile(parameters.profile_reduction_n, profile);
end
# aerosol_optics[iBand][iAer]
aerosol_optics = [Array{vSmartMOM.Scattering.AerosolOptics}(undef, (n_aer)) for i=1:n_bands];
FT2 = isnothing(parameters.scattering_params) ? parameters.float_type : typeof(parameters.scattering_params.rt_aerosols[1].τ_ref)
#FT2 =  parameters.float_type
# τ_aer[iBand][iAer,iZ]
τ_aer = [zeros(FT2, n_aer, length(profile.p_full)) for i=1:n_bands];
kλ_ref = zeros(FT2, n_aer)  # Extinction cross-section at reference wavelength for each aerosol type
kλ     = zeros(FT2, n_aer, n_bands)  # Extinction cross-section at band wavelengths for each aerosol type
ϖλ     = zeros(FT2, n_aer, n_bands)  # Single scattering albedo at band wavelengths for each aerosol type
# Loop over aerosol type
for i_aer=1:n_aer
    # Get curr_aerosol
    c_aero = parameters.scattering_params.rt_aerosols[i_aer]
    curr_aerosol = c_aero.aerosol
    # Create Aerosol size distribution for each aerosol species
    size_distribution = curr_aerosol.size_distribution
    # Create a univariate aerosol distribution
    mie_aerosol = vSmartMOM.Scattering.Aerosol(size_distribution, curr_aerosol.nᵣ, curr_aerosol.nᵢ)
    #@show typeof(curr_aerosol.nᵣ)
    #mie_aerosol = make_mie_aerosol(size_distribution, curr_aerosol.nᵣ, curr_aerosol.nᵢ, parameters.scattering_params.r_max, parameters.scattering_params.nquad_radius) #Suniti: why is the refractive index needed here?
    # Create the aerosol extinction cross-section at the reference wavelength:
    mie_model      = vSmartMOM.Scattering.make_mie_model(parameters.scattering_params.decomp_type,
                                    mie_aerosol,
                                    parameters.scattering_params.λ_ref,
                                    parameters.polarization_type,
                                    truncation_type,
                                    parameters.scattering_params.r_max,
                                    parameters.scattering_params.nquad_radius)
    k_ref          = vSmartMOM.Scattering.compute_ref_aerosol_extinction(mie_model, parameters.float_type)
    kλ_ref[i_aer] = k_ref
    #parameters.scattering_params.rt_aerosols[i_aer].p₀, parameters.scattering_params.rt_aerosols[i_aer].σp
    # Loop over bands
    for i_band=1:n_bands
        # i'th spectral band (convert from cm⁻¹ to μm)
        curr_band_λ = 1e4 ./ parameters.spec_bands[i_band]
        # Create the aerosols:
        mie_model      = vSmartMOM.Scattering.make_mie_model(parameters.scattering_params.decomp_type,
                                        mie_aerosol,
                                        (maximum(curr_band_λ)+minimum(curr_band_λ))/2,
                                        parameters.polarization_type,
                                        truncation_type,
                                        parameters.scattering_params.r_max,
                                        parameters.scattering_params.nquad_radius)
        # Compute raw (not truncated) aerosol optical properties (not needed in RT eventually)
        # @show FT2
        aerosol_optics_raw = vSmartMOM.Scattering.compute_aerosol_optical_properties(mie_model, FT2);
        kλ[i_aer, i_band] = aerosol_optics_raw.k
        ϖλ[i_aer, i_band] = aerosol_optics_raw.ω̃
        # Compute truncated aerosol optical properties (phase function and fᵗ), consistent with Ltrunc:
        #@show i_aer, i_band
        ###aerosol_optics[i_band][i_aer] = vSmartMOM.Scattering.truncate_phase(truncation_type, aerosol_optics_raw; reportFit=false)
        # Compute nAer aerosol optical thickness profiles
        #τ_aer[i_band][i_aer,:] =
        #    parameters.scattering_params.rt_aerosols[i_aer].τ_ref *
        #    (aerosol_optics[i_band][i_aer].k/k_ref) *
        #    vSmartMOM.CoreRT.getAerosolLayerOptProp(1, c_aero.p₀, c_aero.σp, profile.p_half)
        #@show vSmartMOM.CoreRT.getAerosolLayerOptProp(1, c_aero.p₀, c_aero.σp, profile.p_half)
        #@show aerosol_optics[i_band][i_aer].k, k_ref
        #@show τ_aer[i_band][i_aer,:]
        #@show parameters.scattering_params.rt_aerosols[i_aer].τ_ref
    end
end
wl = [770., 2036., 2120., 2204., 2288., 2372.]
τ_all = cat(0.3*ones(5,1), 0.3*kλ./kλ_ref, dims=2)
lstyle = [:dash, :dash, :solid, :solid, :solid]
using Makie
using CairoMakie
# Make a figure
f = Figure(resolution=(800, 600))
ax = Axis(f[1, 1], title="Aerosol within the Carbon-I spectral band (τ at 770 nm = 0.3)",
          xlabel="Wavelength [nm]", ylabel="AOD, τ")
#aer_type = ["water", "sulphate", "sea salt", "mineral dust", "organic"]
aer_type2 = [L"$r_g =0.08\,\mu$m, $\sigma_g=1.6$", L"$r_g =0.13\,\mu$m, $\sigma_g=1.6$", L"$r_g =0.08\,\mu$m, $\sigma_g=1.3$", L"$r_g =0.13\,\mu$m, $\sigma_g=1.3$", L"$r_g =0.20\,\mu$m, $\sigma_g=1.3$"]
# Plot each row
for i=1:size(τ_all, 1)
    lines!(ax, wl[2:end], τ_all[i, 2:end], linestyle=lstyle[i], label=aer_type2[i], linewidth=2)
end
axislegend(ax, position=:rt)
f
save("tau_vs_wl_2.png", f)
f = Figure(resolution=(800, 600))
ax = Axis(f[1, 1], title="Aerosol within the Carbon-I spectral band",
          xlabel="Wavelength [nm]", ylabel="SSA, ϖ")
#aer_type = ["water", "sulphate", "sea salt", "mineral dust", "organic"]
aer_type2 = [L"$r_g =0.08\,\mu$m, $\sigma_g=1.6$", L"$r_g =0.13\,\mu$m, $\sigma_g=1.6$", L"$r_g =0.08\,\mu$m, $\sigma_g=1.3$", L"$r_g =0.13\,\mu$m, $\sigma_g=1.3$", L"$r_g =0.20\,\mu$m, $\sigma_g=1.3$"]
# Plot each row
for i=1:size(τ_all, 1)
    lines!(ax, wl[2:end], ϖλ[i, :], linestyle=lstyle[i], label=aer_type2[i], linewidth=2)
end
axislegend(ax, position=:lb)
f