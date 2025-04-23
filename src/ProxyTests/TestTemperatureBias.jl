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

R_conv_carbonI_dict = Dict{Tuple{Float64, Float64,Float64}, Vector{Float64}}()



parameters = parameters_from_yaml("/home/cfranken/code/gitHub/CarbonI/src/yaml/carbon-i.yaml")

# Let's set the scene for a tropical surface albedo
wl_ = 1e7./parameters.spec_bands[1]
ν = parameters.spec_bands[1];
Δν = mean(diff(ν));
gridi = ν[1]:Δν:ν[end]+10eps()

tProfile = deepcopy(parameters.T)

#

# 10 Degree biases in Boundary layer:
T_shifts = -10:10

n = length(T_shifts)
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
for (iT,dT) in enumerate(T_shifts)
    parameters.T .= tProfile
    # Change lowest 3 layers:
    parameters.T[32:34] .+= dT 
    # Run code:
    model = model_from_parameters(parameters);
    a = rt_run(model)
    rad_inter = CubicSplineInterpolation(gridi, a[1][1,1,:]);
    y = cM*rad_inter(1e7./wl);
    
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
    @show dT, R_ch4./R_n2o, R_ch4
    ratios_ch4[iT]   = R_ch4
    ratios_ch4_3[iT] = R_ch4_3
    ratios_ch4_1[iT] = R_ch4_1
    ratios_ch4_4[iT] = R_ch4_4
    ratios_n2o[iT]   = R_n2o
    mw1_h2o[iT]   = R_h2o_1
    mw2_h2o[iT]   = R_h2o_2
    mw3_h2o[iT]   = R_h2o_3
    mw4_h2o[iT]   = R_h2o_4
    ratios_n2o_w1[iT]   = R_n2o_1
end






function plotTBias()
    f = Figure(resolution=(500,400), title="Impact of wrong prior temperature profile", fontsize=16)

    # Primary axis (left Y)
    ax1 = Axis(f[1, 1];
        xlabel="Boundary Layer temperature bias (K)",
        ylabel="Retrieved XCH₄ using N₂O Proxy method (ppb)",xticks=-10:5:10,
        yminorgridvisible=true,
        #backgroundcolor=:transparent
    )

    

    # plot on ax1
    lines!(ax1,
    T_shifts,
    ratios_ch4./ratios_n2o*2003.5,
        label = "Retrieved XCH₄",
        color = :black,
        linewidth = 3,
    )

    lines!(ax1,
    [-10, 10],
    [2000, 2000],
        label = "True XCH₄",
        color = CarbonI_colors[1],
        linewidth = 3,
    )

    
    
    # separate legends so they don’t overlap
    axislegend(ax1, position = :rc)
  
    CairoMakie.xlims!(ax1, -10, 10)
    
   # xlims!(ax2, 2000, 2400)
    return f
end

f = plotTBias()
save("plots/ProxyPaper/T-Bias-Plots.pdf", f)