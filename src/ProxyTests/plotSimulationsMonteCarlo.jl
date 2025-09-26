using CairoMakie
set_theme!(theme_ggplot2())

using Random
using JLD2, Statistics, vSmartMOM, InstrumentOperator, CarbonI, ImageFiltering, Interpolations, Statistics, CUDA, Plots, Distributions
include("./src/ProxyTests/setupTestFits_generic.jl")
include("./src/ProxyTests/createMicroWindows_cbe.jl")
include("./src/Plots/CI_colors.jl")

Δwl = 0.005
wl = 2035:Δwl:2385

# Define an instrument:
cM, wl_ci = CarbonI.create_carbonI_conv_matrix_cbe(wl)


parameters = parameters_from_yaml("/home/cfranken/code/gitHub/CarbonI/src/yaml/carbon-i.yaml")
aod = 0.1
p_aero = 900.0
albedo = 0.1
sza = 30.0

parameters.scattering_params.rt_aerosols[1].profile = Normal(p_aero, 50.0)
parameters.sza = sza
model = model_from_parameters(parameters);
aod_profile_normalized = model.τ_aer[1][1,:] / sum(model.τ_aer[1][1,:])
#model.τ_rayl[1] .= 0.0
model.τ_aer[1][1,:] = aod * aod_profile_normalized
model.params.brdf[1] = vSmartMOM.CoreRT.LambertianSurfaceScalar{Float64}(albedo)
a = rt_run(model)
aod_band = sum(model.τ_aer[1]); 
ν = parameters.spec_bands[1];
Δν = mean(diff(ν));
gridi = ν[1]:Δν:ν[end]+10eps()
R = a[1][1,1,:];
rad_inter = CubicSplineInterpolation(gridi, R);
R_conv_carbonI = cM*rad_inter(1e7./wl);

# Monte Carlo runs:
n2o_s = []
for i=1:10000
    RR = R_conv_carbonI .+ randn(length(R_conv_carbonI))*0.0001
    x2, yy2, S2, A2, K2 = ii_mw2(RR[indLR2]);
    push!(n2o_s,h_column2["n2o"]' * x2)
end
x2, yy2, S2, A2, K2 = ii_mw2(R_conv_carbonI[indLR2]);

sigma_pred = sqrt(h_column2["n2o"]' * S2 * h_column2["n2o"])
MC_error = std(n2o_s)

# using CairoMakie
using Statistics
using Distributions

# 1) Normalize errors by expected sigma → z should be ~ N(0,1) if predictions are calibrated
    z = (n2o_s .- mean(n2o_s)) ./ sigma_pred

    # 2) Fit a Normal to z
    μ̂   = mean(z)
    σ̂   = std(z, corrected=true)
    fitN = Normal(μ̂, σ̂)

fig = Figure(resolution=(700,400))
ax = Axis(fig[1,1], xlabel="z = MC errors / σ_predicted", ylabel="density", title="Monte Carlo precision validation for N₂O (N=5000)")
bins = -5:0.25:5
hist!(ax, z, bins=bins, normalization=:pdf)
xs = range(-5, 5; length=400)
lines!(ax, xs, pdf.(Normal(0,1), xs),color=:red, label="Ideal N(0,1)")
lines!(ax, xs, pdf.(fitN, xs),color=:orange,  label="Fitted N(μ̂,σ̂)")
axislegend(ax)
CairoMakie.xlims!(-4, 4)
fig

save("plots/ProxyFitsMC.pdf", fig)
save("plots/ProxyFitsMC.png", fig)