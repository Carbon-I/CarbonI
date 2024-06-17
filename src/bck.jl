heatmap(x, bins[2:end], binsDry1' ./sum(binsDry1, dims=2), colormap=:OrRd, clims=(0,0.2), yscale=:log10)
plot(x, percentiles_dry[:,5],linewidth=3,alpha=1,color=:black, label="median", yscale=:log10, ylims=(1e-4,1))
#plot!(x, cf_mean,linewidth=3,alpha=1, color=:red, label="mean")

plot!(x, percentiles_dry[:,4],linewidth=2,alpha=0.8,color=:black, label="40/60 %ile")
plot!(x, percentiles_dry[:,6],linewidth=2,alpha=0.8,color=:black,label=nothing)
plot!(x, percentiles_dry[:,3],linewidth=2,alpha=0.5,color=:black, label="30/70 %ile")
plot!(x, percentiles_dry[:,7],linewidth=2,alpha=0.5,color=:black,label=nothing)
plot!(x, percentiles_dry[:,2],linewidth=1,alpha=0.3,color=:black, label="20/80 %ile")
plot!(x, percentiles_dry[:,8],linewidth=1,alpha=0.3,color=:black,label=nothing)
plot!(x, percentiles_dry[:,1],linewidth=1,alpha=0.1,color=:black, label="10/90 %ile")
plot!(x, percentiles_dry[:,9],linewidth=1,alpha=0.1,color=:black,label=nothing)
xlims!(0,5000)

plot!(x, percentiles_wet[:,5],linewidth=3,alpha=1,color=:blue, label="median")
#plot!(x, cf_mean,linewidth=3,alpha=1, color=:red, label="mean")

plot!(x, percentiles_wet[:,4],linewidth=2,alpha=0.8,color=:blue, label="40/60 %ile")
plot!(x, percentiles_wet[:,6],linewidth=2,alpha=0.8,color=:blue,label=nothing)
plot!(x, percentiles_wet[:,3],linewidth=2,alpha=0.5,color=:blue, label="30/70 %ile")
plot!(x, percentiles_wet[:,7],linewidth=2,alpha=0.5,color=:blue,label=nothing)
plot!(x, percentiles_wet[:,2],linewidth=1,alpha=0.3,color=:blue, label="20/80 %ile")
plot!(x, percentiles_wet[:,8],linewidth=1,alpha=0.3,color=:blue,label=nothing)
plot!(x, percentiles_wet[:,1],linewidth=1,alpha=0.1,color=:blue, label="10/90 %ile")
plot!(x, percentiles_wet[:,9],linewidth=1,alpha=0.1,color=:blue,label=nothing)


plot(fit(Histogram, filter(is_in_dry,datasets[6]).cloud_free_fraction1, bins), seriestype=:steps, lw=3, lc=:red, label="200m Dry Season", xscale=:log10,xlims=(0.9e-4, 0))
plot!(fit(Histogram, filter(is_in_wet,datasets[6]).cloud_free_fraction1, bins), seriestype=:steps, lw=3, lc=:blue, label="200m Wet Season", xscale=:log10,xlims=(0.9e-4, 0))

plot(fit(Histogram, filter(is_in_dry,datasets[15]).cloud_free_fraction1, bins), seriestype=:steps, lw=3, lc=:red, label="2000m Dry Season", xscale=:log10,xlims=(0.9e-4, 0))
plot!(fit(Histogram, filter(is_in_wet,datasets[15]).cloud_free_fraction1, bins), seriestype=:steps, lw=3, lc=:blue, label="2000m Wet Season", xscale=:log10,xlims=(0.9e-4, 0))


filter(is_in_jul_to_sep, df)


plot(x, percentiles[:,5],linewidth=3,alpha=1,color=:black,ylims=(0.0002,1), label="median")
plot!(x, cf_mean,linewidth=3,alpha=1, color=:red, label="mean")

plot!(x, percentiles[:,4],linewidth=2,alpha=0.8,color=:black, label="40/60 %ile")
plot!(x, percentiles[:,6],linewidth=2,alpha=0.8,color=:black,label=nothing)
plot!(x, percentiles[:,3],linewidth=2,alpha=0.5,color=:black, label="30/70 %ile")
plot!(x, percentiles[:,7],linewidth=2,alpha=0.5,color=:black,label=nothing)
plot!(x, percentiles[:,2],linewidth=1,alpha=0.3,color=:black, label="20/80 %ile")
plot!(x, percentiles[:,8],linewidth=1,alpha=0.3,color=:black,label=nothing)
plot!(x, percentiles[:,1],linewidth=1,alpha=0.1,color=:black, label="10/90 %ile")
plot!(x, percentiles[:,9],linewidth=1,alpha=0.1,color=:black,label=nothing)
xlims!(0,5000)
