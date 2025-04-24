using JLD2, Statistics, vSmartMOM, InstrumentOperator, CarbonI, ImageFiltering, Interpolations, Statistics, CUDA, Distributions
using Plots
#device!(1)
Δwl = 0.005
wl = 2035:Δwl:2385

# Define our standard fitting windows:
include("./src/ProxyTests/setupTestFits_generic.jl")
include("./src/ProxyTests/createMicroWindows_cbe.jl")
include("./src/Plots/CI_colors.jl")
# Define an instrument:
cM, wl_ci = CarbonI.create_carbonI_conv_matrix_cbe(wl)




parameters = parameters_from_yaml("/home/cfranken/code/gitHub/CarbonI/src/yaml/carbon-i-splineAlbedo.yaml")

R_conv_carbonI_dict = Dict{Tuple{Float64, Float64,Float64}, Vector{Float64}}()
aods = [0.002, 0.01, 0.02,0.03, 0.04, 0.05, 0.075, 0.1, 0.125,  0.15, 0.175, 0.2]
#aods = [0.002,  0.05,   0.15]
#aods = [0.075]
alb_scalings = 0.05:0.05:0.45
#alb_scaling = [1.0]
flipSigns = [1.0]
#flipSigns = [1.0,   -1.0]
p_aero = 800.0

# Let's set the scene for a tropical surface albedo
wl_ = 1e7./parameters.spec_bands[1]
scenario = CarbonI.stressing_scenario()
parameters.scattering_params.rt_aerosols[1].profile = Normal(p_aero, 50.0)
#parameters.sza = sza
model = model_from_parameters(parameters);
aod_profile_normalized = model.τ_aer[1][1,:] / sum(model.τ_aer[1][1,:])
# End of Basic setup ########################

# Define runs:
n_iter   = 1
ratios_ch4   = Array{Float64}(undef, (length(aods),length(alb_scalings), length(flipSigns), n_iter ));
ratios_ch4_2   = Array{Float64}(undef, (length(aods),length(alb_scalings), length(flipSigns), n_iter ));
ratios_ch4_3 = Array{Float64}(undef, (length(aods),length(alb_scalings), length(flipSigns), n_iter ));
ratios_n2o   = Array{Float64}(undef, (length(aods),length(alb_scalings), length(flipSigns), n_iter ));
ratios_n2o_w1   = Array{Float64}(undef, (length(aods),length(alb_scalings), length(flipSigns), n_iter ));
mw2_h2o   = Array{Float64}(undef, (length(aods),length(alb_scalings), length(flipSigns), n_iter ));
mw3_h2o   = Array{Float64}(undef, (length(aods),length(alb_scalings), length(flipSigns), n_iter ));
mw1_h2o   = Array{Float64}(undef, (length(aods),length(alb_scalings), length(flipSigns), n_iter ));
mw4_h2o   = Array{Float64}(undef, (length(aods),length(alb_scalings), length(flipSigns), n_iter ));

# Loop over AODs and albedo scalings:
for (iAOD,aod) in enumerate(aods)
    for (iAlb,alb_scaling) in enumerate(alb_scalings)
        for (iSign, flipSign) in enumerate(flipSigns)
            println("Iter for AOD, albedo scaling, sign: ")
            @show aod, alb_scaling, flipSign
            model.τ_aer[1][1,:] = aod * aod_profile_normalized
            
            #model.τ_rayl[1] .=1e-30;
            center_fit = mean(wl_ci[indLR2])
            p =  Polynomials.fit(wl_ci[indLR2].-center_fit, scenario.surface_albedo(wl_ci[indLR2]),2)
            gridi = 1980:30:2400
            spl = CubicSplineInterpolation(gridi, ones(length(gridi))*alb_scaling)
            #p[1] = p[1] * flipSign
            #p[2] = p[2] * flipSign
            #grid = 2000:20:2400
            #spl = CubicSplineInterpolation(2000:20:2400, p.(2000-center_fit:20:2400-center_fit)*alb_scaling)
            model.params.brdf[1] = vSmartMOM.CoreRT.LambertianSurfaceSpline(spl, wl_);
            a = rt_run(model)
            ν = parameters.spec_bands[1];
            Δν = mean(diff(ν));
            gridi = ν[1]:Δν:ν[end]+10eps()
            R = a[1][1,1,:];
            rad_inter = CubicSplineInterpolation(gridi, R);
            R_conv_carbonI = cM*rad_inter(1e7./wl);
            key = (aod, alb_scaling, flipSign); 
            R_conv_carbonI_dict[key] = R_conv_carbonI
            
        end
    end
end

@save "simulated_rads_VariableAlbedoBaseline.jld2" R_conv_carbonI_dict 

# x, yy, S, A, K = ii_mw1(R_conv_carbonI[indLR1]);
# Loop over AODs and albedo scalings:
#  @load "simulated_rads_VariableAlbedoTropics.jld2" R_conv_carbonI_dict
for (iAOD,aod) in enumerate(aods)
    for (iAlb,alb_scaling) in enumerate(alb_scalings)
        for (iSign, flipSign) in enumerate(flipSigns[1])
            key = (aod, alb_scaling, flipSign)
            y1 = R_conv_carbonI_dict[key][indLR1];
            y2 = R_conv_carbonI_dict[key][indLR2];
            y3 = R_conv_carbonI_dict[key][indLR3];
            y4 = R_conv_carbonI_dict[key][indLR4];

            x1, yy1, S1, A1, K1 = ii_mw1(y1);
            x2, yy2, S2, A2, K2 = ii_mw2(y2);
            x3, yy3, S3, A3, K3 = ii_mw3(y3);
            x4, yy4, S4, A4, K4 = ii_mw4(y4);
            #y_ = R_conv_carbonI[indLR2]
            R_ch4 = (x2' * h_column2["ch4"])./(xa2' * h_column2["ch4"])
            R_ch4_3 = (x3' * h_column3["ch4"])./(xa3' * h_column3["ch4"])
            R_ch4_2 = (x4' * h_column4["ch4"])./(xa4' * h_column4["ch4"])
            R_n2o = (x2' * h_column2["n2o"])./(xa2' * h_column2["n2o"])
            R_n2o_1 = (x1' * h_column1["n2o"])./(xa1' * h_column1["n2o"])
            R_h2o_1 = (x1' * h_column1["h2o"])./(xa1' * h_column1["h2o"])
            R_h2o_2 = (x2' * h_column2["h2o"])./(xa2' * h_column2["h2o"])
            R_h2o_3 = (x3' * h_column3["h2o"])./(xa3' * h_column3["h2o"])
            R_h2o_4 = (x4' * h_column4["h2o"])./(xa4' * h_column4["h2o"])
            @show key, R_ch4./R_n2o
            ratios_ch4[iAOD,iAlb,iSign,1]   = R_ch4
            ratios_ch4_3[iAOD,iAlb,iSign,1] = R_ch4_3
            ratios_ch4_2[iAOD,iAlb,iSign,1] = R_ch4_2
            ratios_n2o[iAOD,iAlb,iSign,1]   = R_n2o
            mw1_h2o[iAOD,iAlb,iSign,1]   = R_h2o_1
            mw2_h2o[iAOD,iAlb,iSign,1]   = R_h2o_2
            mw3_h2o[iAOD,iAlb,iSign,1]   = R_h2o_3
            mw4_h2o[iAOD,iAlb,iSign,1]   = R_h2o_4
            ratios_n2o_w1[iAOD,iAlb,iSign,1]   = R_n2o_1
        end
    end
end

key = (0.15, 0.35, 1.0)
y1 = R_conv_carbonI_dict[key][indLR1];
y2 = R_conv_carbonI_dict[key][indLR2];
y3 = R_conv_carbonI_dict[key][indLR3];
y4 = R_conv_carbonI_dict[key][indLR4];

x1, yy1, S1, A1, K1 = ii_mw1(y1);
x2, yy2, S2, A2, K2 = ii_mw2(y2);
x3, yy3, S3, A3, K3 = ii_mw3(y3);
x4, yy4, S4, A4, K4 = ii_mw4(y4);


function plotFits()
    f = Figure(resolution=(700,500), title="Spectral Fits", fontsize=15)

    # Primary axis (left Y)
    ax1 = Axis(f[1, 1]; xlabel="Wavelength (μm)", ylabel="Reflected radiance",title="Noise free retrieval performance", yminorgridvisible=true)

    # plot on ax1
    lines!(ax1,wl_ci,R_conv_carbonI_dict[key],label = "Simulated Spectrum",color = :black,linewidth = 2,)
    lines!(ax1,wl_ci[indLR1],yy1,label = "Fit MW1",color = CarbonI_colors[1],linewidth = 2,)
    lines!(ax1,wl_ci[indLR4],yy4,label = "Fit MW2",color = CarbonI_colors[11],linewidth = 2,)
    lines!(ax1,wl_ci[indLR2],yy2,label = "Fit MW3",color = CarbonI_colors[7],linewidth = 2,)
    lines!(ax1,wl_ci[indLR3],yy3,label = "Fit MW4",color = CarbonI_colors[10],linewidth = 2,)
    

    lines!(ax1,wl_ci[indLR1],10*(y1 .- yy1),color = CarbonI_colors[1],linewidth = 2,)
    lines!(ax1,wl_ci[indLR2],10*(y2 .- yy2),color = CarbonI_colors[7],linewidth = 2,)
    lines!(ax1,wl_ci[indLR3],10*(y3 .- yy3),color = CarbonI_colors[10],linewidth = 2,)
    lines!(ax1,wl_ci[indLR4],10*(y4 .- yy4),color = CarbonI_colors[11],linewidth = 2,)
    text!(ax1, "Residuals 10x(measured - fitted)", position = (2200, 0.003), color = :black)
    # separate legends so they don’t overlap
    axislegend(ax1, position = :lb, orientation = :horizontal)
    
    CairoMakie.xlims!(ax1, 2030, 2380)
    CairoMakie.ylims!(ax1, -0.02, 0.6)
    
   # xlims!(ax2, 2000, 2400)
    return f
end

f = plotFits()
save("plots/ProxyPaper/SpectralFitsBaseline.pdf", f)

function plotResults()
    f = Figure(resolution=(800,550), title="Fit Results", fontsize=17)
    ranger = (0.95,1.02)
    # Primary axis (left Y)
    ax1 = Axis(f[1, 1]; xlabel="AOD", ylabel="Reflectance at 2.1µm", yminorgridvisible=true, title="Ω(CH₄)", xticks=0:0.1:0.2)
    ax2 = Axis(f[1, 2]; xlabel="AOD",  yminorgridvisible=true, title="Ω(N₂O)", xticks=0:0.1:0.2)
    ax3 = Axis(f[1, 4]; xlabel="AOD",  yminorgridvisible=true, title="Ω(CH₄)/Ω(N₂O)", xticks=0:0.1:0.2)
    hideydecorations!(ax2, grid=false)
    hideydecorations!(ax3, grid=false)
    ch4 = CairoMakie.heatmap!(ax1,aods,  alb_scalings,ratios_ch4[:,:,1,1], colormap = :viridis, colorbar = false, colorrange=ranger)
    n2o = CairoMakie.heatmap!(ax2,aods, alb_scalings,ratios_n2o[:,:,1,1], colormap = :viridis, colorbar = false, colorrange=ranger)
    proxyRatio = CairoMakie.heatmap!(ax3,aods,  alb_scalings,ratios_ch4[:,:,1,1]./ratios_n2o[:,:,1,1], colormap = :viridis, colorbar = false)
    Colorbar(f[1, 3], ch4)
    Colorbar(f[1, 5], proxyRatio)
    ax4 = Axis(f[2, 1:2]; xlabel="Ω(CH₄, MW3)/Ω(CH₄ MW4)", ylabel="Ω(CH₄, MW3)/Ω(N₂O MW3)", yminorgridvisible=true)
    CairoMakie.scatter!(ax4, (ratios_ch4[:,:,1,1]'./ratios_ch4_3[:,:,1,1]')[:], (ratios_ch4[:,:,1,1]'./ratios_n2o[:,:,1,1]')[:], color = :black)
    ax5 = Axis(f[2, 3:5]; xlabel="Ω(CH₄, MW3)/Ω(CH₄ MW4)", ylabel="Ω(H₂O, MW3)/Ω(H₂O MW4)", yminorgridvisible=true)
    CairoMakie.scatter!(ax5, (ratios_ch4[:,:,1,1]'./ratios_ch4_3[:,:,1,1]')[:], (mw2_h2o[:,:,1,1]'./mw3_h2o[:,:,1,1]')[:], color = :black)
    return f
end
plotResults()
f = plotResults()
save("plots/ProxyPaper/ProxyBaseline2D.pdf", f)