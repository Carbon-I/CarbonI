using vSmartMOM, InstrumentOperator, CarbonI, ImageFiltering, Interpolations, Statistics, CUDA

# Define wavelength grid for forward model:
Δwl = 0.005
wl = 2035:Δwl:2385

# Define an instrument:
cM_CarbonI, wl_carbonI = CarbonI.create_carbonI_conv_matrix(wl)
cM_EMIT, wl_emit    = CarbonI.create_emit_conv_matrix(wl)

# Load parametes and create model
parameters = parameters_from_yaml("src/yaml/carbon-i.yaml")
model = model_from_parameters(parameters);

# Edit parameters here:
model.params.brdf[1] = vSmartMOM.CoreRT.LambertianSurfaceScalar{Float64}(0.05)
ref = sum(model.τ_aer[1][1,:]);
model.τ_aer[1][1,:] .*= 0.1/ref

# run model
a = rt_run(model)

ν = parameters.spec_bands[1];
gridi = LinRange(ν[1], ν[end], length(ν))

# Get reflected radiance:
R = a[1][1,1,:];

# Interpolate to the grid:
rad_inter = CubicSplineInterpolation(gridi, R);
R_wl = rad_inter(1e7./wl);
# Convolve with instrument:
R_conv_carbonI = cM_CarbonI * R_wl
R_conv_emit =    cM_EMIT * R_wl

plot(wl, R_wl, label="HighRes", alpha=0.2)
plot!(wl_carbonI, R_conv_carbonI, label="Carbon-I")
plot!(wl_emit, R_conv_emit, label="EMIT")



