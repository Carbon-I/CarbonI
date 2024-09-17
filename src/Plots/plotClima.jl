using GeoMakie, CairoMakie, NCDatasets 
using Printf

cl = Dataset("/home/wyujie/DATASERVER/model/CLIMA/LAND/simulations/a14_gm2_wd1_2019_1X_1M.hs.epar.kdtd.oldphi.aggpp.int.nc")

cmap = cgrad(:viridis,  15, categorical = true)

lat = cl["lat"][:];
lon = cl["lon"][:];
mGPP = cl["mGPP"][:,:,:];
mSIF740 = cl["mSIF740"][:,:,:];
mSIF757 = cl["mSIF757"][:,:,:];
mPPAR  = cl["mPPAR"][:,:,:];
mΦ_F = cl["mΦ_F"][:,:,:];
mΦ_P = cl["mΦ_P"][:,:,:];
mBETA = cl["mBETA"][:,:,:];

# Plot ClimaYield
i = 6
for i = 1:12
    fig = Figure(resolution=(700,700))
    ga = GeoAxis(title="SIF Yield, $i 2019", 
        fig[1, 1]; # any cell of the figure's layout
        #dest = "+proj=moll",# the CRS in which you want to plot
    )

    ga2 = GeoAxis(title="PSII Yield, $i 2019",
        fig[2, 1]; # any cell of the figure's layout
        #dest = "+proj=moll",# the CRS in which you want to plot
    )
    #sp = GeoMakie.heatmap!(ga, lon, lat, (sum(N, dims=1)[1,:,:]); colorrange=(10^1,10^5), colormap=cmap)
    sp = GeoMakie.heatmap!(ga, lon, lat, mΦ_F[:,:,i]*100;  colormap=cmap, colorrange = (0.8,1.6))
    sp2 = GeoMakie.heatmap!(ga2, lon, lat, mΦ_P[:,:,i] ;  colorrange = (0.1,0.85))

    lines!(ga, GeoMakie.coastlines(), color=:white, alpha=0.5) 
    lines!(ga2, GeoMakie.coastlines(), color=:white, alpha=0.5)
    cb = Colorbar(fig[1, 2], sp; label = "SIF Yield (%)", height = Relative(0.9))
    cb = Colorbar(fig[2, 2], sp2; label = "PSII Yield", height = Relative(0.9))
    rowgap!(fig.layout,1)
    colgap!(fig.layout,1)
    fig
    filename = @sprintf("plots/ClimaYields_%02d.pdf", i)
    CairoMakie.save(filename, fig)
end

for i = 1:12
    fig = Figure(resolution=(700,700))
    ga = GeoAxis(title="SIF @ 757nm, Month= $i", 
        fig[1, 1]; # any cell of the figure's layout
        #dest = "+proj=moll",# the CRS in which you want to plot
    )

    ga2 = GeoAxis(title="GPP Month=$i",
        fig[2, 1]; # any cell of the figure's layout
        #dest = "+proj=moll",# the CRS in which you want to plot
    )
    #sp = GeoMakie.heatmap!(ga, lon, lat, (sum(N, dims=1)[1,:,:]); colorrange=(10^1,10^5), colormap=cmap)
    sp = GeoMakie.heatmap!(ga, lon, lat, mSIF757[:,:,i];  colormap=cmap, colorrange = (0.,0.5))
    sp2 = GeoMakie.heatmap!(ga2, lon, lat, mGPP[:,:,i] ;  colorrange = (0,13))

    lines!(ga, GeoMakie.coastlines(), color=:white, alpha=0.5) 
    lines!(ga2, GeoMakie.coastlines(), color=:white, alpha=0.5)
    cb = Colorbar(fig[1, 2], sp; label = "SIF (mW m⁻² sr⁻¹ nm⁻¹)", height = Relative(0.9))
    cb = Colorbar(fig[2, 2], sp2; label = "GPP (µmol/m²/s)", height = Relative(0.9))
    rowgap!(fig.layout,1)
    colgap!(fig.layout,1)
    fig
    filename = @sprintf("plots/ClimaSIFGPP_%02d.pdf", i)
    CairoMakie.save(filename, fig)
end



for i = 1:12
    a = mGPP[:,:,i]
    a[a .>= 0.1] .= 1.0
    a[a .< 0.1]  .= NaN
    
    fig = Figure(resolution=(700,350))
    ga = GeoAxis(title="Electrons per CO₂, $i 2019", 
        fig[1, 1]; # any cell of the figure's layout
        #dest = "+proj=moll",# the CRS in which you want to plot
    )

    
    #sp = GeoMakie.heatmap!(ga, lon, lat, (sum(N, dims=1)[1,:,:]); colorrange=(10^1,10^5), colormap=cmap)
    sp = GeoMakie.heatmap!(ga, lon, lat, a .* (1e6*(mPPAR[:,:,i] .* mΦ_P[:,:,i])./mGPP[:,:,i])  ;  colormap=cmap, colorrange = (12,50))
    #sp2 = GeoMakie.heatmap!(ga2, lon, lat, mΦ_P[:,:,i] ;  colorrange = (0.1,0.85))

    lines!(ga, GeoMakie.coastlines(), color=:white, alpha=0.5) 
    #lines!(ga2, GeoMakie.coastlines(), color=:white, alpha=0.5)
    cb = Colorbar(fig[1, 2], sp; label = "e-/CO₂ ", height = Relative(0.9))
    #cb = Colorbar(fig[2, 2], sp2; label = "PSII Yield", height = Relative(0.9))
    #rowgap!(fig.layout,1)
    #colgap!(fig.layout,1)
    fig
    filename = @sprintf("plots/ClimaEtoCO2_%02d.pdf", i)
    CairoMakie.save(filename, fig)
end

for i = 1:12
    fig = Figure(resolution=(700,700))
    ga = GeoAxis(title="SIF Yield, $i 2019", 
        fig[1, 1]; # any cell of the figure's layout
        #dest = "+proj=moll",# the CRS in which you want to plot
    )

    ga2 = GeoAxis(title="PSII Yield, $i 2019",
        fig[2, 1]; # any cell of the figure's layout
        #dest = "+proj=moll",# the CRS in which you want to plot
    )
    #sp = GeoMakie.heatmap!(ga, lon, lat, (sum(N, dims=1)[1,:,:]); colorrange=(10^1,10^5), colormap=cmap)
    sp = GeoMakie.heatmap!(ga, lon, lat, mΦ_F[:,:,i]*100;  colormap=cmap, colorrange = (0.8,1.6))
    sp2 = GeoMakie.heatmap!(ga2, lon, lat, mΦ_P[:,:,i] ;  colorrange = (0.1,0.85))

    lines!(ga, GeoMakie.coastlines(), color=:white, alpha=0.5) 
    lines!(ga2, GeoMakie.coastlines(), color=:white, alpha=0.5)
    cb = Colorbar(fig[1, 2], sp; label = "SIF Yield (%)", height = Relative(0.9))
    cb = Colorbar(fig[2, 2], sp2; label = "PSII Yield", height = Relative(0.9))
    rowgap!(fig.layout,1)
    colgap!(fig.layout,1)
    fig
    filename = @sprintf("plots/ClimaYields_%02d.pdf", i)
    CairoMakie.save(filename, fig)
end

for i=1:12
    # Create a figure and axis
    fig = Figure()
    ax = Axis(fig[1, 1], title = "SIF vs GPP, Month=$i", xlabel = "SIF (mW m⁻² sr⁻¹ nm⁻¹)", ylabel = "GPP (μmol CO₂ m⁻² s⁻¹)")

    # Create the scatter plot with color and colormap
    scatterplot = scatter!(ax, mSIF757[:,:,i][:], mGPP[:,:,i][:], color = ((1e6*(mPPAR[:,:,i][:] .* mΦ_P[:,:,i][:])./mGPP[:,:,i][:])), markersize=3, colormap = :viridis, colorrange=(12,50))

    # Set axis limits
    # colorrange=(0.2,0.6)
    xlims!(0, 0.5)
    ylims!(0, 13)


    # Add a colorbar
    Colorbar(fig[1, 2], scatterplot, label = "e⁻/CO₂", height = Relative(0.9))

    # Display the plot
    fig
    filename = @sprintf("plots/ClimaScatter_%02d.pdf", i)
    CairoMakie.save(filename, fig)
end