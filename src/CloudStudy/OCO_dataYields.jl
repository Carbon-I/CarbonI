using NCDatasets, GeoMakie, CairoMakie, MakieThemes
ds = Dataset("/net/fluo/data3/pooled/OCO2/LtCO2/B11014Ar/reprocessed/monthly_1degree/OCO_XCO2-v11.1deg_regrid_GlobalAll.nc");
lat = ds["lat"][:]
lon = ds["lon"][:]
N   = ds["n"][:]
Makie.set_theme!(ggthemr(:earth))

cmap = cgrad(:viridis,  15, categorical = true)
fig = Figure(resolution=(800,400))
ga = GeoAxis(title="Entire OCO-2 data record (2014-2023)",
    fig[1, 1]; # any cell of the figure's layout
    #dest = "+proj=moll",# the CRS in which you want to plot
    coastlines = true # plot coastlines from Natural Earth, as a reference.
)
#sp = GeoMakie.heatmap!(ga, lon, lat, (sum(N, dims=1)[1,:,:]); colorrange=(10^1,10^5), colormap=cmap)
sp = GeoMakie.heatmap!(ga, lon, lat, log10.(sum(N, dims=1)[1,:,:]); colorrange=(1,5), colormap=cmap)
cb = Colorbar(fig[1, 2], sp; label = "Number of Observations", height = Relative(0.85))
cb.axis.attributes[:scale][] = log10
cb.axis.attributes[:limits][] = exp10.(cb.axis.attributes[:limits][])
fig
CairoMakie.save("CarbonI_Clouds_all.pdf", fig)

labels = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
for i in 1:12
    fig = Figure(resolution=(800,400))
    ga = GeoAxis(title=labels[i],
    fig[1, 1]; # any cell of the figure's layout
    #dest = "+proj=moll",# the CRS in which you want to plot
    coastlines = true # plot coastlines from Natural Earth, as a reference.
    )
    sp = GeoMakie.heatmap!(ga, lon, lat,log10.(sum(N[i+4:12:end,:,:], dims=1)[1,:,:]); colorrange=(1,4), colormap=cmap)
    cb = Colorbar(fig[1, 2], sp; label = "Number of Observations", height = Relative(0.85))
    cb.axis.attributes[:scale][] = log10
    cb.axis.attributes[:limits][] = exp10.(cb.axis.attributes[:limits][])
    #cb.axis.ticklabels[] = ["10¹", "10²", "10³"]
    #cb.axis.ticklabels[][1] = "10¹"
    #cb.axis.ticklabels[][2] = "10²"
    #cb.axis.ticklabels[][3] = "10³"
    #cb.axis.attributes[:scale][] = log10
    #cb.axis.attributes[:limits][] = exp10.(cb.axis.attributes[:limits][])
    fig
    CairoMakie.save("CarbonI_Clouds_"*string(i,pad=2)*".pdf", fig)
end

# convert -density 200 -delay 75 CarbonI_Clouds_??.pdf CloudAnims.gifi
