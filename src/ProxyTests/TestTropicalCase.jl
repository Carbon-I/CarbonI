using JLD2, Statistics, vSmartMOM, InstrumentOperator, CarbonI, ImageFiltering, Interpolations, Statistics, CUDA, Distributions
using Plots, Distributions
#device!(1)
Δwl = 0.005
wl = 2035:Δwl:2385

# Define our standard fitting windows:
include("./src/ProxyTests/setupTestFits_generic.jl")
include("./src/ProxyTests/createMicroWindows_cbe.jl")
include("./src/Plots/CI_colors.jl")
# Define an instrument:
cM, wl_ci = CarbonI.create_carbonI_conv_matrix_cbe(wl)


R_conv_carbonI_dict = Dict{Tuple{Float64, Float64,Float64}, Vector{Float64}}()

parameters = parameters_from_yaml("/home/cfranken/code/gitHub/CarbonI/src/yaml/carbon-i-splineAlbedo.yaml")

aod = 0.002

p_aero = 800.0

# Let's set the scene for a tropical surface albedo
wl_ = 1e7./parameters.spec_bands[1]
scenario = CarbonI.stressing_scenario()
#parameters.sza = sza
model = model_from_parameters(parameters);
aod_profile_normalized = model.τ_aer[1][1,:] / sum(model.τ_aer[1][1,:])
# End of Basic setup ########################
model.τ_aer[1][1,:] = aod * aod_profile_normalized
grid = 1980:20:2400
spl = CubicSplineInterpolation(grid, scenario.surface_albedo(grid))
model.params.brdf[1] = vSmartMOM.CoreRT.LambertianSurfaceSpline(spl, wl_);
a = rt_run(model)
ν = parameters.spec_bands[1];
Δν = mean(diff(ν));
gridi = ν[1]:Δν:ν[end]+10eps()
R = a[1][1,1,:];
rad_inter = CubicSplineInterpolation(gridi, R);
R_conv_carbonI = cM*rad_inter(1e7./wl);
            
            
y1 = R_conv_carbonI[indLR1];
y2 = R_conv_carbonI[indLR2];

x1, yy1, S1, A1, K1 = ii_mw1(y1);
x2, yy2, S2, A2, K2 = ii_mw2(y2);

            #y_ = R_conv_carbonI[indLR2]
R_ch4 = (x2' * h_column2["ch4"])./(xa2' * h_column2["ch4"])
#R_ch4_4 = (x4' * h_column3["ch4"])./(xa3' * h_column3["ch4"])
R_n2o = (x2' * h_column2["n2o"])./(xa2' * h_column2["n2o"])
R_n2o_1 = (x1' * h_column1["n2o"])./(xa1' * h_column1["n2o"])
            
@show R_ch4./R_n2o, R_ch4, R_n2o, R_n2o_1

# 1) Pre‑allocate
n_iter = 1_000
ratios     = Vector{Float64}(undef, n_iter)
R_ch4_arr  = Vector{Float64}(undef, n_iter)
R_n2o_arr  = Vector{Float64}(undef, n_iter)

for i=1:n_iter
    y22 = y2 .+ randn(length(y2))*0.0001
    x2, yy2, S2, A2, K2 = ii_mw2(y22);
    R_ch4 = (x2' * h_column2["ch4"])./(xa2' * h_column2["ch4"])
    #R_ch4_4 = (x4' * h_column3["ch4"])./(xa3' * h_column3["ch4"])
    R_n2o = (x2' * h_column2["n2o"])./(xa2' * h_column2["n2o"])
    #R_n2o_1 = (x1' * h_column1["n2o"])./(xa1' * h_column1["n2o"])
    ratios[i] = R_ch4 / R_n2o
    R_ch4_arr[i] = R_ch4
    R_n2o_arr[i] = R_n2o
    @show R_ch4./R_n2o, R_ch4, R_n2o
end

key = (0.1, 1.0, 1.0)
y1 = R_conv_carbonI_dict[key][indLR1];
y2 = R_conv_carbonI_dict[key][indLR2];
y3 = R_conv_carbonI_dict[key][indLR3];
y4 = R_conv_carbonI_dict[key][indLR4];

x1, yy1, S1, A1, K1 = ii_mw1(y1);
x2, yy2, S2, A2, K2 = ii_mw2(y2);
x3, yy3, S3, A3, K3 = ii_mw3(y3);
x4, yy4, S4, A4, K4 = ii_mw4(y4);


function plotFits()
    f = Figure(resolution=(700,500), title="Spectral Fits", fontsize=16)

    # Primary axis (left Y)
    ax1 = Axis(f[1, 1]; xlabel="Wavelength (μm)", ylabel="Reflected radiance", yminorgridvisible=true)


    

    # plot on ax1
    lines!(ax1,wl_ci,R_conv_carbonI_dict[key],label = "Simulated Spectrum",color = :black,linewidth = 2,)
    lines!(ax1,wl_ci[indLR1],yy1,label = "Fit MW1",color = CarbonI_colors[1],linewidth = 2,)
    lines!(ax1,wl_ci[indLR2],yy2,label = "Fit MW2",color = CarbonI_colors[7],linewidth = 2,)
    lines!(ax1,wl_ci[indLR3],yy3,label = "Fit MW3",color = CarbonI_colors[10],linewidth = 2,)
    lines!(ax1,wl_ci[indLR4],yy4,label = "Fit MW4",color = CarbonI_colors[11],linewidth = 2,)

    lines!(ax1,wl_ci[indLR1],10*(y1 .- yy1),color = CarbonI_colors[1],linewidth = 2,)
    lines!(ax1,wl_ci[indLR2],10*(y2 .- yy2),color = CarbonI_colors[7],linewidth = 2,)
    lines!(ax1,wl_ci[indLR3],10*(y3 .- yy3),color = CarbonI_colors[10],linewidth = 2,)
    lines!(ax1,wl_ci[indLR4],10*(y4 .- yy4),color = CarbonI_colors[11],linewidth = 2,)

    # separate legends so they don’t overlap
    axislegend(ax1, position = :lb, orientation = :horizontal)
    
    CairoMakie.xlims!(ax1, 2030, 2380)
    CairoMakie.ylims!(ax1, -0.015, 0.07)
    
   # xlims!(ax2, 2000, 2400)
    return f
end
