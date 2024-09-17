using JLD2, CairoMakie
@load "mw2_fits_all.jld2"
@load "simulated_rads_all.jld2" R_conv_carbonI_dict
szas = 20.0:20:60.0
p_aeros = 750.0:100:950.0

size(szas,1) * size(unique(albedos),1) * size(p_aeros,1) * size(unique(aods),1)
alb_ = reshape(albedos,(size(szas,1), size(p_aeros,1), size(unique(aods),1),size(unique(albedos),1)))
sorted_keys = sort(collect(keys(R_conv_carbonI_dict)));

paods = [a[2] for a in sorted_keys];
szas = [a[1] for a in sorted_keys];

ind_lowlowAOD = findall(paods .== 850 .&& szas .==40 .&& 0.007 .< aods .< 0.01 .&& albedos .> 0.035 )
ind_lowAOD    = findall(paods .== 850 .&& szas .==40 .&& 0.04 .< aods .< 0.042 .&& albedos .> 0.035 )
ind_highAOD   = findall(paods .== 850 .&& szas .==40 .&& 0.16 .< aods .< 0.2 .&& albedos .> 0.035) 

ind_ = findall(paods .== 850 .&& szas .==40 .&& aods .< 0.3 .&& 0.6 .> albedos .> 0.035 )

fig = Figure()
ax = Axis(fig[1, 1], title = "Retrieved/True Trace Gas Ratio",
xlabel = "Albedo", ylabel = "Retrieved/True Ratio", yautolimitmargin = (0.05, 0.15))
# Plot XN₂O lines
line1 = lines!(ax, albedos[ind_lowlowAOD], n2o_mw2[ind_lowlowAOD], label = "XN₂O; AOD = 0.01")
line2 = lines!(ax, albedos[ind_lowAOD], n2o_mw2[ind_lowAOD], label = "XN₂O; AOD = 0.04")
line3 = lines!(ax, albedos[ind_highAOD], n2o_mw2[ind_highAOD], label = "XN₂O; AOD = 0.18")

# Plot CH₄ lines with the same colors as XN₂O but dashed
#lines!(ax, albedos[ind_lowlowAOD], ch4_mw2[ind_lowlowAOD], label = "XCH₄; AOD = 0.01", color = line1[:color], linestyle = :dash)
#lines!(ax, albedos[ind_lowAOD], ch4_mw2[ind_lowAOD], label = "XCH₄; AOD = 0.04", color = line2[:color], linestyle = :dash)
#lines!(ax, albedos[ind_highAOD], ch4_mw2[ind_highAOD], label = "XCH₄; AOD = 0.18", color = line3[:color], linestyle = :dash)
xlims!(0.03, 0.6)
ylims!(0.96, 1.04)
axislegend(ax, position = :rb)  # `position = :rb` places the legend at the right bottom
save("plots/n2o_mw2_lowlowAOD_1.pdf", fig)
# Add the legend

fig

fig = Figure()
ax = Axis(fig[1, 1], title = "Retrieved/True Trace Gas Ratio",
xlabel = "Albedo", ylabel = "Retrieved/True Ratio", yautolimitmargin = (0.05, 0.15))
# Plot XN₂O lines
line1 = lines!(ax, albedos[ind_lowlowAOD], n2o_mw2[ind_lowlowAOD], label = "XN₂O; AOD = 0.01")
line2 = lines!(ax, albedos[ind_lowAOD], n2o_mw2[ind_lowAOD], label = "XN₂O; AOD = 0.04")
line3 = lines!(ax, albedos[ind_highAOD], n2o_mw2[ind_highAOD], label = "XN₂O; AOD = 0.18")

# Plot CH₄ lines with the same colors as XN₂O but dashed
lines!(ax, albedos[ind_lowlowAOD], ch4_mw2[ind_lowlowAOD], label = "XCH₄; AOD = 0.01", color = line1[:color], linestyle = :dash)
lines!(ax, albedos[ind_lowAOD], ch4_mw2[ind_lowAOD], label = "XCH₄; AOD = 0.04", color = line2[:color], linestyle = :dash)
lines!(ax, albedos[ind_highAOD], ch4_mw2[ind_highAOD], label = "XCH₄; AOD = 0.18", color = line3[:color], linestyle = :dash)
xlims!(0.03, 0.6)
ylims!(0.96, 1.04)
axislegend(ax, position = :rb)  # `position = :rb` places the legend at the right bottom
save("plots/n2o_ch4_mw2_lowlowAOD_1.pdf", fig)
# Add the legend

fig

fig = Figure()
ax = Axis(fig[1, 1], title = "CH₄/N₂O",
xlabel = "Albedo", ylabel = "Retrieved/True Ratio", yautolimitmargin = (0.05, 0.15))
# Plot XN₂O lines
line1 = lines!(ax, albedos[ind_lowlowAOD], ch4_mw2[ind_lowlowAOD]./n2o_mw2[ind_lowlowAOD], label = "XN₂O; AOD = 0.01")
line2 = lines!(ax, albedos[ind_lowAOD], ch4_mw2[ind_lowAOD]./n2o_mw2[ind_lowAOD], label = "XN₂O; AOD = 0.04")
line3 = lines!(ax, albedos[ind_highAOD], ch4_mw2[ind_highAOD]./n2o_mw2[ind_highAOD], label = "XN₂O; AOD = 0.18")

# Plot CH₄ lines with the same colors as XN₂O but dashed
#lines!(ax, albedos[ind_lowlowAOD], ch4_mw2[ind_lowlowAOD], label = "XCH₄; AOD = 0.01", color = line1[:color], linestyle = :dash)
#lines!(ax, albedos[ind_lowAOD], ch4_mw2[ind_lowAOD], label = "XCH₄; AOD = 0.04", color = line2[:color], linestyle = :dash)
#lines!(ax, albedos[ind_highAOD], ch4_mw2[ind_highAOD], label = "XCH₄; AOD = 0.18", color = line3[:color], linestyle = :dash)
xlims!(0.03, 0.6)
ylims!(0.996, 1.0)
axislegend(ax, position = :rt)  # `position = :rb` places the legend at the right bottom
save("plots/n2o_proxy_ratio.pdf", fig)
# Add the legend

fig


using CairoMakie, NCDatasets, Polynomials

oco2 = Dataset("data/oco2_LtCO2_180825_B11100Ar_230602034736s.nc4")

lat = oco2["latitude"][:]
lon = oco2["longitude"][:]
sza = oco2["solar_zenith_angle"][:] 
xco2_s = oco2.group["Preprocessors"]["xco2_strong_idp"][:]
xco2_w = oco2.group["Preprocessors"]["xco2_weak_idp"][:]
xco2_raw = oco2.group["Retrieval"]["xco2_raw"][:]
xco2 = oco2["xco2"][:]
albedo_s = oco2.group["Retrieval"]["albedo_sco2"][:]
albedo_w = oco2.group["Retrieval"]["albedo_wco2"][:]

ind = findall(52.5 .< lat .< 53 .&& 38.5 .< lon .< 40.25)
ind2 = findall(52.5 .< lat .< 53 .&& 38.5 .< lon .< 40.25 .&& xco2 .< 406)

# Polynomial fit
p = fit(convert.(Float32,albedo_s[ind2]), convert.(Float32,xco2_s[ind2]), 2)

xco2_s_corr = xco2_s[ind] .- p.(convert.(Float32,albedo_s[ind]))


albs = 0.05:0.01:0.3
fig1 = Figure()
ax1 = Axis(fig1[1, 1], title = "Albedo vs XCO₂",
           xlabel = "Albedo", ylabel = "XCO₂ (ppm)")

scatter!(ax1, albedo_s[ind], xco2_s[ind], label = "Strong CO₂ Non-Scattering", color = :blue)
lines!(ax1, albs, p.(convert.(Float32,albs)), label = "Polynomial fit", color = :red)

axislegend(ax1, position = :rb)
display(fig1)

fig2 = Figure()
ax2 = Axis(fig2[1, 1], title = "XCO₂ vs Albedo-Corrected Non-Scattering XCO₂",
           xlabel = "Full Physics XCO₂ - 405 ppm", ylabel = "Albedo-Corrected Non-Scattering XCO₂")

scatter!(ax2, xco2_raw[ind] .- 405, xco2_s_corr, label = "Raw FP", color = :green)
scatter!(ax2, xco2[ind] .- 405, xco2_s_corr, label = "Bias Corrected FP", color = :purple)

axislegend(ax2, position = :rb)
display(fig2)
save("plots/oco2_plume1.pdf", fig1)
save("plots/oco2_plume2.pdf", fig2)

# Create a figure and axis
fig3 = Figure()
ax3 = Axis(fig3[1, 1], title = "Map of XCO₂",
           xlabel = "Longitude", ylabel = "Latitude")

# Scatter plot with color mapping based on xco2 values
scatterplot = scatter!(ax3, lon[ind], lat[ind], color = convert.(Float64,xco2[ind]), colormap = :viridis, 
                       label = "XCO₂ FP Bias-Corrected", markersize = 10,colorrange = (403, 410))

# Add a colorbar
cb = Colorbar(fig3[1, 2], scatterplot, label = "XCO₂ (ppm)")

# Add a legend
axislegend(ax3, position = :rt)

# Display the figure
display(fig3)

# Create a figure and axis
fig4 = Figure()
ax4 = Axis(fig4[1, 1], title = "Map of XCO₂",
           xlabel = "Longitude", ylabel = "Latitude")

# Scatter plot with color mapping based on xco2 values
scatterplot = scatter!(ax4, lon[ind], lat[ind], color = convert.(Float64,xco2_s_corr), colormap = :viridis, 
                       label = "XCO₂ Non-Scattering Bias-Corrected", markersize = 10, colorrange = (-2.5, 4.5))

# Add a colorbar
cb = Colorbar(fig4[1, 2], scatterplot, label = "XCO₂ (ppm)")

# Add a legend
axislegend(ax4, position = :rt)

# Display the figure
display(fig4)
save("plots/oco2_plumeMap1.pdf", fig4)
save("plots/oco2_plumeMap2.pdf", fig3)