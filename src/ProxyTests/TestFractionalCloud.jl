using JLD2, Statistics, vSmartMOM, InstrumentOperator, CarbonI, ImageFiltering, Interpolations, Statistics, CUDA, Distributions
using DelimitedFiles, Plots
#device!(1)
Δwl = 0.005
wl = 2035:Δwl:2385

# Define our standard fitting windows:
include("./src/ProxyTests/setupTestFits_generic.jl")
include("./src/ProxyTests/createMicroWindows_cbe.jl")
include("./src/Plots/CI_colors.jl")
# Define an instrument:
cM, wl_ci = CarbonI.create_carbonI_conv_matrix_cbe(wl)

scenario = CarbonI.stressing_scenario()



grid = 2000:10:2400
spl = CubicSplineInterpolation(grid, scenario.surface_albedo(grid))


parameters = parameters_from_yaml("/home/cfranken/code/gitHub/CarbonI/src/yaml/carbon-i-cloud.yaml")
model_cloud = model_from_parameters(parameters);
model_cloud.τ_aer[1][1,:] .= 0.0
# Put cloud in one layer only 800-850hPa
model_cloud.τ_aer[1][1,31] = 20.0
# Let's set the scene for a tropical surface albedo
wl_ = 1e7./parameters.spec_bands[1]
ν = parameters.spec_bands[1];
Δν = mean(diff(ν));
gridi = ν[1]:Δν:ν[end]+10eps()

# Run with cloud
a_cloud = rt_run(model_cloud)
rad_inter = CubicSplineInterpolation(gridi, a_cloud[1][1,1,:]);
R_conv_cloud = cM*rad_inter(1e7./wl);

model_cloud.τ_aer[1][1,:] .= 0.0
# Run without cloud
a_no_cloud = rt_run(model_cloud)
rad_inter = CubicSplineInterpolation(gridi, a_no_cloud[1][1,1,:]);
R_conv_NoCloud = cM*rad_inter(1e7./wl);

cloud_fractions = 0:0.01:1

n = length(cloud_fractions)
# Define runs:

ratios_ch4   = Array{Float64}(undef, (n));
ratios_ch4_3 = Array{Float64}(undef, (n));
ratios_ch4_1 = Array{Float64}(undef, (n));
ratios_ch4_4 = Array{Float64}(undef, (n));
ratios_n2o   = Array{Float64}(undef, (n));
ratios_n2o_w1   = Array{Float64}(undef, (n));
mw2_h2o   = Array{Float64}(undef, (n));
mw3_h2o   = Array{Float64}(undef, (n));
mw1_h2o   = Array{Float64}(undef, (n));
mw4_h2o   = Array{Float64}(undef, (n));

# Loop over AODs and albedo scalings:
for (icF,cF) in enumerate(cloud_fractions)
    # IPA independent Pixel approximation here:
    y = (1-cF)*R_conv_NoCloud + cF*R_conv_cloud 

    y1 = y[indLR1];
    y2 = y[indLR2];
    y3 = y[indLR3];
    y4 = y[indLR4]
    x1, yy1, S1, A1, K1 = ii_mw1(y1);
    x2, yy2, S2, A2, K2 = ii_mw2(y2);
    x3, yy3, S3, A3, K3 = ii_mw3(y3);
    x4, yy4, S4, A4, K4 = ii_mw4(y4);
    #y_ = R_conv_carbonI[indLR2]
    R_ch4 = (x2' * h_column2["ch4"])./(xa2' * h_column2["ch4"])
    R_ch4_1 = (x1' * h_column1["ch4"])./(xa1' * h_column1["ch4"])
    R_ch4_4 = (x4' * h_column4["ch4"])./(xa4' * h_column4["ch4"])
    R_ch4_3 = (x3' * h_column3["ch4"])./(xa3' * h_column3["ch4"])
    R_n2o = (x2' * h_column2["n2o"])./(xa2' * h_column2["n2o"])
    R_n2o_1 = (x1' * h_column1["n2o"])./(xa1' * h_column1["n2o"])
    R_h2o_1 = (x1' * h_column1["h2o"])./(xa1' * h_column1["h2o"])
    R_h2o_2 = (x2' * h_column2["h2o"])./(xa2' * h_column2["h2o"])
    R_h2o_3 = (x3' * h_column3["h2o"])./(xa3' * h_column3["h2o"])
    R_h2o_4 = (x4' * h_column4["h2o"])./(xa4' * h_column4["h2o"])
    @show cF, R_ch4./R_n2o, R_ch4
    ratios_ch4[icF]   = R_ch4
    ratios_ch4_3[icF] = R_ch4_3
    ratios_ch4_1[icF] = R_ch4_1
    ratios_ch4_4[icF] = R_ch4_4
    ratios_n2o[icF]   = R_n2o
    mw1_h2o[icF]   = R_h2o_1
    mw2_h2o[icF]   = R_h2o_2
    mw3_h2o[icF]   = R_h2o_3
    mw4_h2o[icF]   = R_h2o_4
    ratios_n2o_w1[icF]   = R_n2o_1
end




function plotCloudFractional()
    f = Figure(resolution=(500,400), title="Impact of fractional cloud cover on boundary layer enhancements", fontsize=16)

    # Primary axis (left Y)
    ax1 = Axis(f[1, 1];
        xlabel="Cloud Fraction (%)",
        ylabel="Retrieved XCH₄ using N₂O Proxy method (ppb)",title="Impact of thick shallow cumulus cloud",
        yminorgridvisible=true,xticks=0:10:100,
        #backgroundcolor=:transparent
    )

    

    # plot on ax1
    lines!(ax1,
    cloud_fractions*100,
    ratios_ch4./ratios_n2o*2000,
        label = "Retrieved XCH₄",
        color = :black,
        linewidth = 3,
    )

    lines!(ax1,
    [0, 100],
    [2046.5, 2046.5],
        label = "True total column XCH₄",
        color = CarbonI_colors[1],
        linewidth = 3,
    )

    lines!(ax1,
    [0, 100],
    [2000, 2000],
        label = "True above-cloud XCH₄",
        color = CarbonI_colors[10],
        linewidth = 3,
    )
    
    # separate legends so they don’t overlap
    axislegend(ax1, position = :lb)
  
    CairoMakie.xlims!(ax1, 0, 100)
    CairoMakie.ylims!(ax1, 1950, 2070)
    
   # xlims!(ax2, 2000, 2400)
    return f
end

f = plotCloudFractional()
save("plots/ProxyPaper/FractionalClouds.pdf", f)