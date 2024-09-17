using CarbonI, GeoMakie, CairoMakie

geos, aero = CarbonI.loadGeos("src/GeosChem/GeosChem.yaml");

dp = CarbonI.computeColumnAveragingOperator(geos);
tropo_lev = geos.data["tropopause_level"];
lat = geos.data["lat"];
lon = geos.data["lon"].+2.5;
xn2o = CarbonI.getColumnAverage(geos.data["N2O"], dp);
xch4 = CarbonI.getColumnAverage(geos.data["CH4"], dp);
xch4_trop = CarbonI.getTroposphericColumnAverage(geos.data["CH4"], dp, floor.(tropo_lev));
xn2o_trop = CarbonI.getTroposphericColumnAverage(geos.data["N2O"], dp, floor.(tropo_lev));

c2h6 = geos.data["C2H6"];
xc2h6 = CarbonI.getColumnAverage(c2h6, dp)
cmap = cgrad(:viridis,  15, categorical = true)

# Plot CH4:
fig = Figure(resolution=(700,700))
ga = GeoAxis(title="Tropospheric Column Average, June 2019",
    fig[1, 1]; # any cell of the figure's layout
    #dest = "+proj=moll",# the CRS in which you want to plot
)

ga2 = GeoAxis(title="Total Column Average, June 2019",
    fig[2, 1]; # any cell of the figure's layout
    #dest = "+proj=moll",# the CRS in which you want to plot
)
#sp = GeoMakie.heatmap!(ga, lon, lat, (sum(N, dims=1)[1,:,:]); colorrange=(10^1,10^5), colormap=cmap)
sp = GeoMakie.contourf!(ga, lon, lat, xch4_trop*1e9;  colormap=cmap)
sp2 = GeoMakie.contourf!(ga2, lon, lat, xch4*1e9;  colormap=cmap)

lines!(ga, GeoMakie.coastlines(), color=:white, alpha=0.5) 
lines!(ga2, GeoMakie.coastlines(), color=:white, alpha=0.5)
cb = Colorbar(fig[1, 2], sp; label = "XCH₄ (ppb) ", height = Relative(0.9))
cb = Colorbar(fig[2, 2], sp2; label = "XCH₄ (ppb) ", height = Relative(0.9))
rowgap!(fig.layout,1)
colgap!(fig.layout,1)
fig
CairoMakie.save("plots/CarbonI_XCH4.pdf", fig)
# Plot N2O:
fig = Figure(resolution=(700,700))
ga = GeoAxis(title="Tropospheric Column Average, June 2019",
    fig[1, 1]; # any cell of the figure's layout
    #dest = "+proj=moll",# the CRS in which you want to plot
)

ga2 = GeoAxis(title="Total Column Average, June 2019",
    fig[2, 1]; # any cell of the figure's layout
    #dest = "+proj=moll",# the CRS in which you want to plot
)
#sp = GeoMakie.heatmap!(ga, lon, lat, (sum(N, dims=1)[1,:,:]); colorrange=(10^1,10^5), colormap=cmap)
sp = GeoMakie.contourf!(ga, lon, lat, xn2o_trop*1e9;  colormap=cmap)
sp2 = GeoMakie.contourf!(ga2, lon, lat, xn2o*1e9;  colormap=cmap)

lines!(ga, GeoMakie.coastlines(), color=:white, alpha=0.5) 
lines!(ga2, GeoMakie.coastlines(), color=:white, alpha=0.5)
cb = Colorbar(fig[1, 2], sp; label = "N₂O (ppb) ", height = Relative(0.9))
cb = Colorbar(fig[2, 2], sp2; label = "N₂O (ppb) ", height = Relative(0.9))
rowgap!(fig.layout,1)
colgap!(fig.layout,1)
fig
CairoMakie.save("plots/CarbonI_XN2O.pdf", fig)
# Plot N2O:
fig = Figure(resolution=(700,700))
ga = GeoAxis(title="CH₄ Tropospheric/Total Column Average, June 2019",
    fig[1, 1]; # any cell of the figure's layout
    #dest = "+proj=moll",# the CRS in which you want to plot
)

ga2 = GeoAxis(title="N₂O Tropospheric/Total Column Average, June 2019",
    fig[2, 1]; # any cell of the figure's layout
    #dest = "+proj=moll",# the CRS in which you want to plot
)
#sp = GeoMakie.heatmap!(ga, lon, lat, (sum(N, dims=1)[1,:,:]); colorrange=(10^1,10^5), colormap=cmap)
sp = GeoMakie.contourf!(ga, lon, lat, xch4_trop./xch4;  colormap=cmap)
sp2 = GeoMakie.contourf!(ga2, lon, lat, xn2o_trop./xn2o;  colormap=cmap)

lines!(ga, GeoMakie.coastlines(), color=:white, alpha=0.5) 
lines!(ga2, GeoMakie.coastlines(), color=:white, alpha=0.5)
cb = Colorbar(fig[1, 2], sp; label = "Trop/Total ratio ", height = Relative(0.8))
cb = Colorbar(fig[2, 2], sp2; label = "Trop/Total ratio ", height = Relative(0.8))
rowgap!(fig.layout,1)
colgap!(fig.layout,1)
fig
CairoMakie.save("plots/CarbonI_tropColumnRatio.pdf", fig)

# Plot N2O:
fig = Figure(resolution=(700,350))
ga = GeoAxis(title="Tropopause height, June 2019",
    fig[1, 1]; # any cell of the figure's layout
    #dest = "+proj=moll",# the CRS in which you want to plot
)


#sp = GeoMakie.heatmap!(ga, lon, lat, (sum(N, dims=1)[1,:,:]); colorrange=(10^1,10^5), colormap=cmap)
sp = GeoMakie.contourf!(ga, lon, lat, geos.data["tropopause_height"][:,:,1];  colormap=cmap)


lines!(ga, GeoMakie.coastlines(), color=:white, alpha=0.5) 
#lines!(ga2, GeoMakie.coastlines(), color=:white, alpha=0.5)
cb = Colorbar(fig[1, 2], sp; label = "Tropopause height (km) ", height = Relative(0.8))
#cb = Colorbar(fig[2, 2], sp2; label = "Trop/Total ratio ", height = Relative(0.8))
rowgap!(fig.layout,1)
colgap!(fig.layout,1)
fig
CairoMakie.save("plots/CarbonI_tropopauseH.pdf", fig)
