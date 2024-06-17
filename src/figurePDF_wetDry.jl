cmap = cgrad(:YlGnBu,  15, categorical = true)

# Create a Figure
fig = Figure(resolution=(800,400))

# Create an axis with a logarithmic y-scale
ax = Axis(fig[1, 1], yscale = log10,title = "Amazon Dry Season")

data = (binsDry1' ./sum(binsDry1, dims=2)')'
data_ = log10.(data .+ 1e-10)
hm = CairoMakie.heatmap!(ax, x, [0.8e-4; bins[2:end]], data_, colormap=cmap,ylims=(1e-4,1),xlims = (0, 2000), colorrange=(-1.5,-0))# Add a colorbar
#cb = Colorbar(fig[1, 2], hm, label = "Probability (log10)", labelpadding = 10, ticklabelpad = 5, width = 20)
#cb.axis.attributes[:scale][] = log10
#cb.axis.attributes[:limits][] = exp10.(cb.axis.attributes[:limits][])


# Overlay the line plot
iVal = findall(percentiles_dry[:,5] .>1e-4)
lines!(ax, x[iVal], percentiles_dry[iVal,5],linewidth=3,alpha=1,color=:black, label="median",ylims=(1e-4,1),xlims = (0, 2000))
lines!(ax, x[iVal], percentiles_dry[iVal,6],linewidth=2,alpha=0.7,color=:black, label="median",ylims=(1e-4,1),xlims = (0, 2000))
lines!(ax, x[iVal], percentiles_dry[iVal,7],linewidth=2,alpha=0.5,color=:black, label="median",ylims=(1e-4,1),xlims = (0, 2000))
lines!(ax, x[iVal], percentiles_dry[iVal,8],linewidth=1,alpha=0.3,color=:black, label="median",ylims=(1e-4,1),xlims = (0, 2000))
lines!(ax, x[iVal], percentiles_dry[iVal,9],linewidth=1,alpha=0.1,color=:black, label="median",ylims=(1e-4,1),xlims = (0, 2000))
iVal = findall(percentiles_dry[:,4] .>1e-4)
lines!(ax, x[iVal], percentiles_dry[iVal,4],linewidth=2,alpha=0.7,color=:black, label="median",ylims=(1e-4,1),xlims = (0, 2000))
iVal = findall(percentiles_dry[:,3] .>1e-4)
lines!(ax, x[iVal], percentiles_dry[iVal,3],linewidth=2,alpha=0.5,color=:black, label="median",ylims=(1e-4,1),xlims = (0, 2000))
iVal = findall(percentiles_dry[:,2] .>1e-4)
lines!(ax, x[iVal], percentiles_dry[iVal,2],linewidth=1,alpha=0.3,color=:black, label="median",ylims=(1e-4,1),xlims = (0, 2000))
iVal = findall(percentiles_dry[:,1] .>1e-4)
lines!(ax, x[iVal], percentiles_dry[iVal,1],linewidth=1,alpha=0.1,color=:black, label="median",ylims=(1e-4,1),xlims = (0, 2000))


# Show the figure
#fig
#CairoMakie.save("Cloud_DrySeason_S2.pdf", fig)

# Create a Figure
#fig = Figure()

# Create an axis with a logarithmic y-scale
ax2 = Axis(fig[1, 2], yscale = log10, title = "Amazon Wet Season")

# Plot the heatmap
data = (binsWet1' ./sum(binsWet1, dims=2)')'
data_ = log10.(data .+ 1e-10)
hm = CairoMakie.heatmap!(ax2, x, [0.8e-4; bins[2:end]], data_,colormap=cmap, ylims=(0.8e-4,1),xlims = (0, 2500), colorrange=(-1.5,0))
# Add a colorbar
#cb.axis.attributes[:scale][] = log10
cb = Colorbar(fig[1, 3], hm,  label = "Probability (log10)", labelpadding = 10, ticklabelpad = 5, width = 20)
#cb.axis.attributes[:scale][] = log10
#cb.axis.attributes[:limits][] = exp10.(cb.axis.attributes[:limits][])
# Overlay the line plot
iVal = findall(percentiles_wet[:,5] .>1e-4)
lines!(ax2, x[iVal], percentiles_wet[iVal,5],linewidth=3,alpha=1,color=:black, label="median",ylims=(1e-4,1),xlims = (0, 2500))
lines!(ax2, x[iVal], percentiles_wet[iVal,6],linewidth=2,alpha=0.7,color=:black, label="median",ylims=(1e-4,1),xlims = (0, 2500))
lines!(ax2, x[iVal], percentiles_wet[iVal,7],linewidth=2,alpha=0.5,color=:black, label="median",ylims=(1e-4,1),xlims = (0, 2500))
lines!(ax2, x[iVal], percentiles_wet[iVal,8],linewidth=1,alpha=0.3,color=:black, label="median",ylims=(1e-4,1),xlims = (0, 2500))
lines!(ax2, x[iVal], percentiles_wet[iVal,9],linewidth=1,alpha=0.1,color=:black, label="median",ylims=(1e-4,1),xlims = (0, 2500))
iVal = findall(percentiles_wet[:,4] .>1e-4)
lines!(ax2, x[iVal], percentiles_wet[iVal,4],linewidth=2,alpha=0.7,color=:black, label="median",ylims=(1e-4,1),xlims = (0, 2500))
iVal = findall(percentiles_wet[:,3] .>1e-4)
lines!(ax2, x[iVal], percentiles_wet[iVal,3],linewidth=2,alpha=0.5,color=:black, label="median",ylims=(1e-4,1),xlims = (0, 2500))
iVal = findall(percentiles_wet[:,2] .>1e-4)
lines!(ax2, x[iVal], percentiles_wet[iVal,2],linewidth=1,alpha=0.3,color=:black, label="median",ylims=(1e-4,1),xlims = (0, 2500))
iVal = findall(percentiles_wet[:,1] .>1e-4)
lines!(ax2, x[iVal], percentiles_wet[iVal,1],linewidth=1,alpha=0.1,color=:black, label="median",ylims=(1e-4,1),xlims = (0, 2500))


# Show the figure
colsize!(fig.layout, 1, Relative(0.47))
colsize!(fig.layout, 2, Relative(0.47))
#colsize!(fig.layout, 3, Relative(0.2))
ax2.xlabel = "Footprint size (m)"
ax.xlabel = "Footprint size (m)"
ax.ylabel = "Fraction of cloud free pixels (%)"
ax2.yticklabelsvisible = false
ax.yticks = ([0.0001,0.001,0.01,0.1,1], ["<0.01","0.1","1","10","100"])



ii = 14
# Add a rotated label at the chosen point
xx = x

text!(ax, "Median", position = (xx[ii],percentiles_dry[ii,5]), rotation = -0.17, color = :red)
text!(ax, "60%ile", position = (xx[ii],percentiles_dry[ii,6]), rotation = -0.16, color = :red)
text!(ax, "70%ile", position = (xx[ii],percentiles_dry[ii,7]), rotation = -0.1, color = :red)
text!(ax, "80%ile", position = (xx[ii],percentiles_dry[ii,8]), rotation = -0.08, color = :red)
text!(ax, "90%ile", position = (xx[ii],percentiles_dry[ii,9]), rotation = -0.025, color = :red)
text!(ax, "10%ile", position = (xx[3],percentiles_dry[3,1]+0.0001), rotation = -1.3, color = :red)
text!(ax, "20%ile", position = (xx[11],percentiles_dry[11,2]), rotation = -0.75, color = :red)
text!(ax, "30%ile", position = (xx[ii],percentiles_dry[ii,3]), rotation = -0.4, color = :red)
text!(ax, "40%ile", position = (xx[ii],percentiles_dry[ii,4]), rotation = -0.25, color = :red)

# Add a rotated label at the chosen point
text!(ax2, "Median", position = (xx[ii],percentiles_wet[ii,5]), rotation = -0.6, color = :red)
text!(ax2, "60%ile", position = (xx[ii],percentiles_wet[ii,6]), rotation = -0.4, color = :red)
text!(ax2, "70%ile", position = (xx[ii],percentiles_wet[ii,7]), rotation = -0.3, color = :red)
text!(ax2, "80%ile", position = (xx[ii],percentiles_wet[ii,8]), rotation = -0.2, color = :red)
text!(ax2, "90%ile", position = (xx[ii],percentiles_wet[ii,9]), rotation = -0.08, color = :red)
#text!(ax2, "10%ile", position = (xx[3],percentiles_wet[3,1]+0.0001), rotation = -1.3, color = :red)
#text!(ax2, "20%ile", position = (xx[11],percentiles_wet[11,2]), rotation = -0.75, color = :red)
text!(ax2, "30%ile", position = (xx[6],percentiles_wet[6,3]), rotation = -0.9, color = :red)
text!(ax2, "40%ile", position = (xx[11],percentiles_wet[11,4]), rotation = -0.7, color = :red)


fig
CairoMakie.save("Cloud_WetDrySeason_S2.pdf", fig)

function compAngle(xx,yy,index)
    x1 = xx[index-1]
    x2 = xx[index]
    y1 = yy[index-1]
    y2 = yy[index]
    return atan((y2-y1)/(x2-x1))
end