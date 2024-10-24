using NCDatasets, GeoMakie, CairoMakie, Statistics, Dates

mod = Dataset("/net/fluo/data1/DATASERVER/satellite/MODIS/MCD43A4.006/reprocessed/0.0833_lat-lon_8d/global/modis_MCD43A43-006.refl.00833deg_regrid.8d.2020.nc")

alb_2 = mod["refl_band7"][:]
lat = mod["lat"][:]
lon = mod["lon"][:]
time = mod["time"][:]
levels = 0.02, 0.04, 0.06, 0.1, 0.2, 0.4, 0.6, 1.0
d = size(alb_2)
alb_avg_summer = zeros(d[2],d[3])
alb_avg_winter = zeros(d[2],d[3])
alb_all = zeros(d[2],d[3])
summer_indices = findall(dt -> month(dt) in (6, 7, 8), time)
winter_indices = findall(dt -> month(dt) in (12, 1, 2), time)

for i=1:d[2]
    for j=1:d[3]
        n = 0;
        for k in summer_indices
            if !ismissing(alb_2[k,i,j]) &&  alb_2[k,i,j] >= 0.025 &&  alb_2[k,i,j] <= 1
                alb_avg_summer[i,j] += alb_2[k,i,j]
                n += 1
            end
        end
        alb_avg_summer[i,j] /= n

        n = 0;
        for k in winter_indices
            if !ismissing(alb_2[k,i,j]) &&  alb_2[k,i,j] >= 0.025 &&  alb_2[k,i,j] <= 1
                alb_avg_winter[i,j] += alb_2[k,i,j]
                n += 1
            end
        end
        alb_avg_winter[i,j] /= n

        n = 0;
        for k in 1:d[1]
            if !ismissing(alb_2[k,i,j]) &&  alb_2[k,i,j] >= 0.025 &&  alb_2[k,i,j] <= 1
                alb_all[i,j] += alb_2[k,i,j]
                n += 1
            end
        end
        alb_all[i,j] /= n
    end
end



clima_alb = readdlm("data/albedo.csv",',', skipstart=1)
modis_limits = [2105.0,2155.0];

Plots.plot(clima_alb[:,1],clima[:,2])
Plots.plot(clima_alb[:,1],clima[:,2])

using CairoMakie
# Create a figure and axis
fig = Figure(resolution=(650,500))

# Create an axis with a logarithmic y-scale
ax = Axis(fig[1, 1], xlabel="Wavelength (nm)",ylabel="Albedo", title = "Tropical Forest Albedo")
CairoMakie.lines!(ax,clima_alb[:,1],clima_alb[:,2]/1.16, linewidth=2, color=:black, label="CliMA land model") 

# Define the rectangle corners
xs = [modis_limits[1], modis_limits[2], modis_limits[2], modis_limits[1]]  # x-coordinates of the rectangle
ys = [0.0, 0.0, 0.1, 0.1]  # y-coordinates of the rectangle

# Add a transparent rectangle
CairoMakie.poly!(ax, xs, ys, color = (:blue, 0.1), label="MODIS Band 7")  # 30% transparency
CairoMakie.xlims!(2030, 2380)
CairoMakie.ylims!(0, 0.1)
# Create a legend and add it to the figure
axislegend(ax, position=:rt)
# Display the plot
fig
save("plots/tropics_albedo.pdf", fig)

#using MakieThemes
#Makie.set_theme!(ggthemr(:earth))

# Define levels and corresponding colors
levels = [0.0, 0.2, 0.5, 0.7, 1.0] # non-linear levels
colors = [:blue, :green, :yellow, :orange, :red] # corresponding colors

# Create a custom colormap
colormap = colorscheme(colors, levels)

scaling = let thr = 2.0, lthr = log(thr)
    function scale(x)
      (x > thr) ? log(x) : ((x < -thr) ? -log(-x) : (x/thr)*lthr)
    end
    function invscale(x)
      (x > lthr) ? exp(x) : ((x < -lthr) ? -exp(-x) : (x/lthr)*thr)
    end
    CairoMakie.ReversibleScale(scale, invscale)
end

cmap = cgrad(:viridis,  10, categorical = true)
levels = [0.0, 2, 5, 7, 10,15, 20, 40, 60, 80,100]
cScale = LinearInterpolation(levels, collect(1:length(levels)), extrapolation_bc=NaN)

fig = Figure(resolution=(600,400))
ga = GeoAxis(title="MODIS Surface Albedo at 2.15µm",
    fig[1, 1]; 
)
#sp = GeoMakie.heatmap!(ga, lon, lat, (sum(N, dims=1)[1,:,:]); colorrange=(10^1,10^5), colormap=cmap)
sp = GeoMakie.heatmap!(ga, lon, lat, cScale.(100*alb_all); colormap=cmap)
cb = Colorbar(fig[1, 2], sp; label = "Albedo (%)", height = Relative(0.85))
cb.ticks[] = (collect(1:length(levels)), string.(levels))
lines!(ga, GeoMakie.coastlines(), alpha=0.3) 
fig
CairoMakie.save("CarbonI_Clouds_all.png", fig)

indi = findall(abs.(lat).<55)
alb_all[:,indi][:]
alb_all_flt = alb_all[alb_all.>0.01]

fig = Figure(resolution=(450,300))
# Create an axis with a logarithmic y-scale
ax = Axis(fig[1, 1], xlabel="Albedo",ylabel="Frequency", title = "Albedo distribution at 2.1µm")
CairoMakie.hist!(ax,alb_all_flt,color = :red,bins = 150,  normalization=:pdf)
lines!(ax, [0.06, 0.06], [0, 13], color = :black, linewidth = 2, linestyle = :dash)
lines!(ax, [median(alb_all_flt), median(alb_all_flt)], [0, 13], color = :orange, linewidth = 2, linestyle = :dash)
lines!(ax, [mean(alb_all_flt), mean(alb_all_flt)], [0, 13], color = :orange, linewidth = 2)

median(alb_all_flt)
mean(alb_all_flt)
CairoMakie.xlims!(0.02, 0.5)
CairoMakie.ylims!(0, 13)
fig
CairoMakie.save("CarbonI_Albedo_dist.pdf", fig)
N₂O