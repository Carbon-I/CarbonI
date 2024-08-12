using NCDatasets, CairoMakie, GeoMakie, ImageFiltering

# Load Edgar Data:
ed_ch4 = Dataset("data/edgar/v8.0_FT2022_GHG_CH4_2022_TOTALS_flx.nc");
ed_co2 = Dataset("data/edgar/v8.0_FT2022_GHG_CO2_2022_TOTALS_flx.nc");
ed_n2o = Dataset("data/edgar/v8.0_FT2022_GHG_N2O_2022_TOTALS_flx.nc");

# Get Lat/lon
lat = ed_ch4["lat"][:];
lon = ed_ch4["lon"][:];

# Get Fluxes:
flux_ch4 = ed_ch4["fluxes"][:];
flux_co2 = ed_co2["fluxes"][:];
flux_n2o = ed_n2o["fluxes"][:];

# Smoothen the fields as extreme values won't plot well in a single 10km pixel
KernelSize = 10;
kernel = ones(10,10); 
kernel ./= sum(kernel)
step = floor(Int, KernelSize/2)
flux_ch4_smooth = imfilter(flux_ch4, kernel);
flux_co2_smooth = imfilter(flux_co2, kernel);
flux_n2o_smooth = imfilter(flux_n2o, kernel);

year_to_seconds = 3.154e+7;

# Example Plots
fig = Figure(size = (1400, 800))
ga = GeoAxis(fig[1,1])
ch4_field = heatmap!(ga, lon[1:step:end],lat[1:step:end],1000*year_to_seconds*flux_ch4_smooth[1:step:end,1:step:end]; colorrange=(0.01, 100.0), colormap = :OrRd_9)
lines!(ga, GeoMakie.coastlines()) # plot coastlines from Natural Earth as a reference
cb1 = Colorbar(fig[1,2], ch4_field; label = "CH₄ emissions (g/m²/year)", height = Relative(0.65))
fig
save("plots/ch4_fluxes_2022.png", fig; px_per_unit=2)

# CO2:
fig = Figure(size = (1400, 800))
ga = GeoAxis(fig[1,1])
co2_field = heatmap!(ga, lon[1:step:end],lat[1:step:end],year_to_seconds*flux_co2_smooth[1:step:end,1:step:end]; colorrange=(0.1, 10.0), colormap = :OrRd_9)
lines!(ga, GeoMakie.coastlines()) # plot coastlines from Natural Earth as a reference
cb1 = Colorbar(fig[1,2], co2_field; label = "CO₂ emissions (kg/m²/year)", height = Relative(0.65))
fig
save("plots/co2_fluxes_2022.png", fig; px_per_unit=2)

# N2O:
fig = Figure(size = (1400, 800))
ga = GeoAxis(fig[1,1])
co2_field = heatmap!(ga, lon[1:step:end],lat[1:step:end],1000*year_to_seconds*flux_n2o_smooth[1:step:end,1:step:end]; colorrange=(0.01, 2.0), colormap = :OrRd_9)
lines!(ga, GeoMakie.coastlines()) # plot coastlines from Natural Earth as a reference
cb1 = Colorbar(fig[1,2], co2_field; label = "N₂O emissions (g/m²/year)", height = Relative(0.65))
fig
save("plots/n2o_fluxes_2022.png", fig; px_per_unit=2)