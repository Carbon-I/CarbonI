# From silumations:
using JLD2, CairoMakie, ColorSchemes, LaTeXStrings, Polynomials

CS = ColorSchemes.seaborn_colorblind
include(joinpath(@__DIR__, "src/Plots", "CI_colorsNew.jl"))
# Load spectra and retrievals
@load "simulated_rads_all.jld2" R_conv_carbonI_dict
@load "mw2_fits_all.jld2"
@load "mw1_fits_all.jld2"

sorted_keys = sort(collect(keys(R_conv_carbonI_dict)));

paods = [a[2] for a in sorted_keys];
szas = [a[1] for a in sorted_keys];

albedos = convert.(Float64, albedos)
aods = convert.(Float64, aods)
n2o_mw1_ = convert.(Float64,n2o_mw1);
co2_mw1_ = convert.(Float64,co2_mw1);
ch4_mw2_ = convert.(Float64,ch4_mw2);
n2o_mw2_ = convert.(Float64,n2o_mw2);

# polynomial fits for N2O vs CO2:
ind_ = findall(albedos .> 0.04 .&& aods .< 0.5)
p = fit(n2o_mw2_[ind_],co2_mw1_[ind_],2)

# Find some pairs with low, middle and high AOD
ind_lowlowAOD = findall(paods .== 850 .&& szas .==40 .&& 0.0083 .< aods .< 0.011 .&& albedos .> 0.035 )
ind_lowAOD = findall(paods .== 850 .&& szas .==40 .&& 0.045 .< aods .< 0.05 .&& albedos .> 0.035 )
ind_highAOD = findall(paods .== 850 .&& szas .==40 .&& 0.16 .< aods .< 0.2 .&& albedos .> 0.035) 


ind_lowAlb = findall(paods .== 950 .&& szas .==40 .&&  0.061 .> albedos .> 0.06  )
ind_highAlb = findall(paods .== 950 .&& szas .==40 .&& 0.25 .> albedos .> 0.24 ) 


f = Figure(resolution=(800,600), title="Aerosol Microphysics", fontsize=16)
lw = 3
ax1 = Axis(f[1,1],xlabel=L"\text{Albedo}",ylabel=L"\text{Retrieval Bias (%)}",  title=L"\text{Scattering induced retrieval bias in \Omega}", xminorgridvisible = true, xminorticks = IntervalsBetween(5))
lines!(ax1, albedos[ind_lowlowAOD], (-1 .+ n2o_mw2_[ind_lowlowAOD])*100, label=L"\text{\Omega_{N_2O} AOD = 0.01}", color=CarbonI_colors[1], linewidth=lw)
lines!(ax1, albedos[ind_lowlowAOD], (-1 .+ ch4_mw2_[ind_lowlowAOD])*100 .+ 0.25, label=L"\text{\Omega_{CH_4} AOD = 0.01}", color=CarbonI_colors[1], linewidth=lw, linestyle=:dash)
lines!(ax1, albedos[ind_lowAOD], (-1 .+ n2o_mw2_[ind_lowAOD])*100, label=L"\text{\Omega_{N_2O} AOD = 0.05}", color=CarbonI_colors[7], linewidth=lw)
lines!(ax1, albedos[ind_lowAOD], (-1 .+ ch4_mw2_[ind_lowAOD])*100 .+ 0.25, label=L"\text{\Omega_{CH_4} AOD = 0.05}", color=CarbonI_colors[7], linewidth=lw, linestyle=:dash)
lines!(ax1, albedos[ind_highAOD], (-1 .+ n2o_mw2_[ind_highAOD])*100, label=L"\text{\Omega_{N_2O} AOD = 0.165}", color=CarbonI_colors[10], linewidth=lw)
lines!(ax1, albedos[ind_highAOD], (-1 .+ ch4_mw2_[ind_highAOD])*100 .+ 0.25, label=L"\text{\Omega_{CH_4} AOD = 0.165}", color=CarbonI_colors[10], linewidth=lw, linestyle=:dash)
CairoMakie.xlims!(ax1, 0.04, 0.7)
xticks_values = [0.1,0.2,0.3,0.4,0.5,0.6,0.7]
xticks_labels = ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7"]

# Apply custom tick labels
ax1.xticks = (xticks_values, xticks_labels)
axislegend(ax1, position=:rb)
CairoMakie.ylims!(ax1, -4.8, 4.8)

ax2 = Axis(f[1,2], xlabel=L"\text{Albedo}", title=L"\text{Scattering induced retrieval bias in proxy ratio}", xminorgridvisible = true, xminorticks = IntervalsBetween(5))
lines!(ax2, albedos[ind_lowlowAOD], (-1 .+ ch4_mw2_[ind_lowlowAOD]./n2o_mw2_[ind_lowlowAOD])*100 .+ 0.25, label=L"\text{\Omega_{CH_4}/\Omega_{N_2O} AOD = 0.01}", color=CarbonI_colors[1], linewidth=lw)
#lines!(ax1, albedos[ind_lowlowAOD], (-1 .+ ch4_mw2_[ind_lowlowAOD])*100, label=L"\text{XCH_4; AOD = 0.01}", color=CS[1], linewidth=2, linestyle=:dash)
lines!(ax2, albedos[ind_lowAOD], (-1 .+ ch4_mw2_[ind_lowAOD]./n2o_mw2_[ind_lowAOD])*100 .+ 0.25, label=L"\text{\Omega_{CH_4}/\Omega_{N_2O} AOD = 0.05}", color=CarbonI_colors[7], linewidth=lw)
#lines!(ax1, albedos[ind_lowAOD], (-1 .+ ch4_mw2_[ind_lowAOD])*100, label=L"\text{XCH_2; AOD = 0.05}", color=CS[2], linewidth=2, linestyle=:dash)
lines!(ax2, albedos[ind_highAOD], (-1 .+ ch4_mw2_[ind_highAOD]./n2o_mw2_[ind_highAOD])*100 .+ 0.25, label=L"\text{\Omega_{CH_4}/\Omega_{N_2O} AOD = 0.165}", color=CarbonI_colors[10], linewidth=lw)
#lines!(ax1, albedos[ind_highAOD], (-1 .+ ch4_mw2_[ind_highAOD])*100, label=L"\text{XCH_4; AOD = 0.165}", color=CS[3], linewidth=2, linestyle=:dash)
CairoMakie.xlims!(ax2, 0.04, 0.7)
ax2.xticks = (xticks_values, xticks_labels)
axislegend(ax2, position=:rt)
#CairoMakie.ylims!(ax2, -0.5, 0.5)

ax3 = Axis(f[2,1],xlabel=L"\text{AOD}",ylabel=L"\text{Retrieval Bias (%)}",  xminorgridvisible = true, xminorticks = IntervalsBetween(5))
#lines!(ax3, albedos[ind_lowlowAOD], (-1 .+ n2o_mw2_[ind_lowlowAOD])*100, label=L"\text{\Omega_{N_2O} AOD = 0.01}", color=CS[1], linewidth=2)
#lines!(ax3, albedos[ind_lowlowAOD], (-1 .+ ch4_mw2_[ind_lowlowAOD])*100, label=L"\text{\Omega_{CH_4} AOD = 0.01}", color=CS[1], linewidth=2, linestyle=:dash)
lines!(ax3, aods[ind_lowAlb], (-1 .+ n2o_mw2_[ind_lowAlb])*100, label=L"\text{\Omega_{N_2O} Albedo = 0.06}", color=CarbonI_colors[7], linewidth=lw)
lines!(ax3, aods[ind_lowAlb], (-1 .+ ch4_mw2_[ind_lowAlb])*100 .+ 0.25, label=L"\text{\Omega_{CH_4} Albedo = 0.06}", color=CarbonI_colors[7], linewidth=lw, linestyle=:dash)
lines!(ax3, aods[ind_highAlb], (-1 .+ n2o_mw2_[ind_highAlb])*100, label=L"\text{\Omega_{N_2O} Albedo = 0.25}", color=CarbonI_colors[10], linewidth=lw)
lines!(ax3, aods[ind_highAlb], (-1 .+ ch4_mw2_[ind_highAlb])*100 .+ 0.25, label=L"\text{\Omega_{CH_4} Albedo = 0.25}", color=CarbonI_colors[10], linewidth=lw, linestyle=:dash)
CairoMakie.xlims!(ax3, 0.01, 0.3)
CairoMakie.ylims!(ax3, -4.8, 4.8)
axislegend(ax3, position=:lb)
#CairoMakie.ylims!(ax3, -4.8, 4.8)

ax4 = Axis(f[2,2], xlabel=L"\text{AOD}", xminorgridvisible = true, xminorticks = IntervalsBetween(5))
#lines!(ax4, albedos[ind_lowlowAOD], (-1 .+ ch4_mw2_[ind_lowlowAOD]./n2o_mw2_[ind_lowlowAOD])*100 .+ 0.25, label=L"\text{\Omega_{CH_4}/\Omega_{N_2O} AOD = 0.01}", color=CS[1], linewidth=2)
#lines!(ax1, albedos[ind_lowlowAOD], (-1 .+ ch4_mw2_[ind_lowlowAOD])*100, label=L"\text{XCH_4; AOD = 0.01}", color=CS[1], linewidth=2, linestyle=:dash)
lines!(ax4, aods[ind_lowAlb], (-1 .+ ch4_mw2_[ind_lowAlb]./n2o_mw2_[ind_lowAlb])*100 .+ 0.25, label=L"\text{\Omega_{CH_4}/\Omega_{N_2O} Albedo = 0.06}", color=CarbonI_colors[7], linewidth=lw)
#lines!(ax1, albedos[ind_lowAOD], (-1 .+ ch4_mw2_[ind_lowAOD])*100, label=L"\text{XCH_2; AOD = 0.05}", color=CS[2], linewidth=2, linestyle=:dash)
lines!(ax4, aods[ind_highAlb], (-1 .+ ch4_mw2_[ind_highAlb]./n2o_mw2_[ind_highAlb])*100 .+ 0.25, label=L"\text{\Omega_{CH_4}/\Omega_{N_2O} Albedo = 0.25}", color=CarbonI_colors[10], linewidth=lw)
#lines!(ax1, albedos[ind_highAOD], (-1 .+ ch4_mw2_[ind_highAOD])*100, label=L"\text{XCH_4; AOD = 0.165}", color=CS[3], linewidth=2, linestyle=:dash)
CairoMakie.xlims!(ax4, 0.01, 0.3)
CairoMakie.ylims!(ax4, -0.25, 0.25)
CairoMakie.ylims!(ax2, -0.25, 0.25)
axislegend(ax4, position=:lt)
#CairoMakie.ylims!(ax4, -0.5, 0.5)


f
save("plots/ProxyRT_errorsCI_3.pdf", f)
save("plots/final/SectionD-D12-ProxyRT_errorsCI_3.eps", f)
save("plots/final/SectionD-D12-ProxyRT_errorsCI_3.png", f)





