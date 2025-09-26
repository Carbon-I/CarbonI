using CairoMakie
using JLD2, Statistics, vSmartMOM, InstrumentOperator, CarbonI 
using ImageFiltering, Interpolations, Statistics, CUDA, Distributions

# Define some colors and themes
set_theme!(CairoMakie.theme_black())
# Define the color using the hex code
my_color = RGB(0xDC / 255, 0xCA / 255, 0xB0 / 255)  # Normalize the hex values to [0, 1]

function create_gaussian_conv_matrix(modelling_wl::StepRangeLen{FT}, instrument_wl, FWHM) where FT
    # Define a Fixed instrument:
    Δwl = modelling_wl.step.hi
    
    # LSF (Gaussian)
    kernf = CarbonI.gaussian_kernel(FWHM, Δwl)
    
    #Get the Instrument-specific box
    lociBox = CarbonI.KernelInstrument(kernf, instrument_wl);

    # Generate final convolution matrix:
    cM = CarbonI.generate_conv_matrix(lociBox,modelling_wl, Δwl)

    return cM, lociBox.ν_out

end

include("./src/Plots/CI_colors.jl")

# Wavelength grid for high resolution simulations:
Δwl = 0.005
wl = 2035:Δwl:2385

# Define Albedo scenario:
# Load stressing scenario (tropical forest)
scenario = CarbonI.stressing_scenario()
cbe_specs = CarbonI.build_instrument("CBE") 

soil_cbe, x_cbe, solarIrr_cbe, σ_matrix_cbe, profile_cbe, h_cbe, Sₐ_cbe = setup_data(scenario, cbe_specs);
refl_cbe   = scenario.surface_albedo(cbe_specs.modelling_wl);

# Define an instrument operator (for Carbon-I and EMIT):
cM, wl_ci = CarbonI.create_carbonI_conv_matrix_cbe(wl)
cM_emit, wl_emit = CarbonI.create_emit_conv_matrix(wl)
cM_oco2, wl_oco2 = create_gaussian_conv_matrix(wl, collect(2035:0.045:2380), 0.10)  # OCO-2 0.1 nm FWHM


# Load default scenario and run the RT model (the fit has the same input, T, H2O and other profiles could be shuffled here to create inconsistencies)
parameters = parameters_from_yaml("/home/cfranken/code/gitHub/CarbonI/src/yaml/carbon-i.yaml")
aod = 0.1
p_aero = 900.0
albedo = 0.1
sza = 30.0

# Update parameters with the variables set above:
parameters.scattering_params.rt_aerosols[1].profile = Normal(p_aero, 50.0)
parameters.sza = sza

# Generate the model from the parameters:
model = model_from_parameters(parameters);
aod_profile_normalized = model.τ_aer[1][1,:] / sum(model.τ_aer[1][1,:])
model.τ_aer[1][1,:] = aod * aod_profile_normalized
# Lambertian albedo, this could be changed to a polynomial or a more complex model as well
model.params.brdf[1] = vSmartMOM.CoreRT.LambertianSurfaceScalar{Float64}(albedo)

# Run the model:
a = rt_run(model)

# Generate the convolved radiances for EMIT and Carbon-I 
ν = parameters.spec_bands[1];
Δν = mean(diff(ν));
gridi = ν[1]:Δν:ν[end]+10eps()
R = a[1][1,1,:];
rad_inter = CubicSplineInterpolation(gridi, R);
##################################################################

F = rad_inter(1e7./wl) .* scenario.surface_albedo(wl);

# Convolve all simulations:
R_conv_carbonI = cM*F;
R_conv_EMIT    = cM_emit*F;
R_conv_OCO2    = cM_oco2*F;


# Plot starts here:
f = Figure(resolution=(600,400),  fontsize=13, backgroundcolor = :black, fontcolor=:white)
ax1 = Axis(f[1,1],yminorgridvisible = true,  yminorticks = IntervalsBetween(5), xlabel="Wavelength (nm)",ylabel="Reflected Radiance",  title="RT Simulations")
lines!(ax1,wl_oco2, 10R_conv_OCO2, color=:orange, linewidth=1.0, alpha=0.5)
CairoMakie.xlims!(ax1, 2040, 2378)
CairoMakie.ylims!(ax1, 0, 0.06)
save("plots/sim_oco2.pdf", f)
#lines!(ax1,wl_ci, R_conv_carbonI, color = :orange,linewidth=3)
#save("plots/Fits_2.pdf", f)


f

lines!(ax1,wl_emit, 10R_conv_EMIT, color = :red,linewidth=3, alpha=0.8)


save("plots/sim_oco2_EMIT.pdf", f)
lines!(ax1,wl_ci, 10R_conv_carbonI, color = :white,linewidth=2, alpha=0.8)
save("plots/sim_oco2_EMIT_carbonI.pdf", f)

f

f = Figure(resolution=(600,400),  fontsize=13, backgroundcolor = :black, fontcolor=:white)
ax1 = Axis(f[1,1],yminorgridvisible = true,  yminorticks = IntervalsBetween(5), xlabel="Wavelength (nm)",ylabel="Reflected Radiance",  title="RT Simulations")
lines!(ax1,wl, rad_inter(1e7./wl), color=my_color, linewidth=0.2, alpha=0.4)
CairoMakie.xlims!(ax1, 2040, 2378)
CairoMakie.ylims!(ax1, 0, 0.1)


save("plots/Fits_3.pdf", f)
lines!(ax1,wl_ci, R_conv_carbonI, color = :orange,linewidth=2.5, alpha=0.75)
save("plots/Fits_4.pdf", f)


f


