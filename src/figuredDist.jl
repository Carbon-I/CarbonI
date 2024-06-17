cmap = cgrad(:YlGnBu,  15, categorical = true)

# Create a Figure
fig = Figure(resolution=(850,500))

# Create an axis with a logarithmic y-scale
ax = Axis(fig[1, 1], yscale = log10,title = "Amazon Dry Season")

data = (binsDry1' ./sum(binsDry1, dims=2)')'
data_ = log10.(data .+ 1e-10)
data_ = data
data_dry = data_[:,1]
hm = CairoMakie.heatmap!(ax, x2, [0.8e-3; bins[2:end]], data_, colormap=cmap,ylims=(1e-3,1),xlims = (0, 2000), colorrange=(0,0.1), alpha=0.7)# Add a colorbar
#cb = Colorbar(fig[1, 2], hm, label = "Probability (log10)", labelpadding = 10, ticklabelpad = 5, width = 20)
#cb.axis.attributes[:scale][] = log10
#cb.axis.attributes[:limits][] = exp10.(cb.axis.attributes[:limits][])
x2 = [x[1:15]; 2275; 2550; 3000]

# Overlay the line plot
iVal = findall(percentiles_dry[:,5] .>1e-4)
lines!(ax, x2[iVal], percentiles_dry[iVal,5],linewidth=3,alpha=1,color=:black, label="median",ylims=(1e-3,1),xlims = (0, 2000))
iVal = findall(percentiles_dry[:,6] .>1e-4)
lines!(ax, x2[iVal], percentiles_dry[iVal,6],linewidth=2,alpha=0.7,color=:black, label="median",ylims=(1e-3,1),xlims = (0, 2000))
iVal = findall(percentiles_dry[:,7] .>1e-4)
lines!(ax, x2[iVal], percentiles_dry[iVal,7],linewidth=2,alpha=0.5,color=:black, label="median",ylims=(1e-3,1),xlims = (0, 2000))
iVal = findall(percentiles_dry[:,8] .>1e-4)
lines!(ax, x2[iVal], percentiles_dry[iVal,8],linewidth=1,alpha=0.3,color=:black, label="median",ylims=(1e-3,1),xlims = (0, 2000))
iVal = findall(percentiles_dry[:,9] .>1e-4)
lines!(ax, x2[iVal], percentiles_dry[iVal,9],linewidth=1,alpha=0.1,color=:black, label="median",ylims=(1e-3,1),xlims = (0, 2000))
iVal = findall(percentiles_dry[:,4] .>1e-4)
lines!(ax, x2[iVal], percentiles_dry[iVal,4],linewidth=2,alpha=0.7,color=:black, label="median",ylims=(1e-3,1),xlims = (0, 2000))
iVal = findall(percentiles_dry[:,3] .>1e-4)
lines!(ax, x2[iVal], percentiles_dry[iVal,3],linewidth=2,alpha=0.5,color=:black, label="median",ylims=(1e-3,1),xlims = (0, 2000))
iVal = findall(percentiles_dry[:,2] .>1e-4)
lines!(ax, x2[iVal], percentiles_dry[iVal,2],linewidth=1,alpha=0.3,color=:black, label="median",ylims=(1e-3,1),xlims = (0, 2000))
iVal = findall(percentiles_dry[:,1] .>1e-4)
lines!(ax, x2[iVal], percentiles_dry[iVal,1],linewidth=1,alpha=0.1,color=:black, label="median",ylims=(1e-3,1),xlims = (0, 2000))

m_line = lines!(ax, x2, cf_mean_dry,linewidth=2,alpha=0.5,color=:red)
# Show the figure
#fig
#CairoMakie.save("Cloud_DrySeason_S2.pdf", fig)

# Create a Figure
#fig = Figure()

# Create an axis with a logarithmic y-scale
ax2 = Axis(fig[1, 2], yscale = log10, title = "Amazon Wet Season")

# Plot the heatmap
data = (binsWet1' ./sum(binsWet1, dims=2)')'
data_ = data#  log10.(data .+ 1e-10)
data_wet = data_[:,1]
hm = CairoMakie.heatmap!(ax2, x2, [0.8e-3; bins[2:end]], data_,colormap=cmap, ylims=(0.8e-4,1),xlims = (0, 2500), colorrange=(0,0.1), alpha=0.7)
# Add a colorbar
#cb.axis.attributes[:scale][] = log10
cb = Colorbar(fig[1, 3], hm,  label = "Probability distribution", labelpadding = 10, ticklabelpad = 5, width = 20)
#cb.axis.attributes[:scale][] = log10
#cb.axis.attributes[:limits][] = exp10.(cb.axis.attributes[:limits][])
# Overlay the line plot
iVal = findall(percentiles_wet[:,5] .>1e-4)
med_line = lines!(ax2, x2[iVal], percentiles_wet[iVal,5],linewidth=3,alpha=1,color=:black, label="median",ylims=(1e-3,1),xlims = (0, 2500))
iVal = findall(percentiles_wet[:,6] .>1e-4)
lines!(ax2, x2[iVal], percentiles_wet[iVal,6],linewidth=2,alpha=0.7,color=:black, label="median",ylims=(1e-3,1),xlims = (0, 2500))
iVal = findall(percentiles_wet[:,7] .>1e-4)
perc1 = lines!(ax2, x2[iVal], percentiles_wet[iVal,7],linewidth=2,alpha=0.5,color=:black, label="median",ylims=(1e-3,1),xlims = (0, 2500))
iVal = findall(percentiles_wet[:,8] .>1e-4)
lines!(ax2, x2[iVal], percentiles_wet[iVal,8],linewidth=1,alpha=0.3,color=:black, label="median",ylims=(1e-3,1),xlims = (0, 2500))
iVal = findall(percentiles_wet[:,9] .>1e-4)
lines!(ax2, x2[iVal], percentiles_wet[iVal,9],linewidth=1,alpha=0.1,color=:black, label="median",ylims=(1e-3,1),xlims = (0, 2500))
iVal = findall(percentiles_wet[:,4] .>1e-4)
lines!(ax2, x2[iVal], percentiles_wet[iVal,4],linewidth=2,alpha=0.7,color=:black, label="median",ylims=(1e-3,1),xlims = (0, 2500))
iVal = findall(percentiles_wet[:,3] .>1e-4)
lines!(ax2, x2[iVal], percentiles_wet[iVal,3],linewidth=2,alpha=0.5,color=:black, label="median",ylims=(1e-3,1),xlims = (0, 2500))
iVal = findall(percentiles_wet[:,2] .>1e-4)
lines!(ax2, x2[iVal], percentiles_wet[iVal,2],linewidth=1,alpha=0.3,color=:black, label="median",ylims=(1e-3,1),xlims = (0, 2500))
iVal = findall(percentiles_wet[:,1] .>1e-4)
lines!(ax2, x2[iVal], percentiles_wet[iVal,1],linewidth=1,alpha=0.1,color=:black, label="median",ylims=(1e-3,1),xlims = (0, 2500))
mean_line = lines!(ax2, x2, cf_mean_wet,linewidth=2,alpha=0.5,color=:red, label="mean")
#axislegend(ax2, [mean_line,med_line,perc1 ], ["Mean", "Median","Percentiles"],  position = :rt,
#    orientation = :horizontal)

Legend(fig[2, 3], [mean_line,med_line,perc1 ], ["Mean", "Median","Percentiles"], framevisible = false)
# Show the figure
colsize!(fig.layout, 1, Relative(0.45))
colsize!(fig.layout, 2, Relative(0.45))
colsize!(fig.layout, 3, Relative(0.1))
#colsize!(fig.layout, 3, Relative(0.2))

ax.ylabel = "Fraction of cloud free pixels (%)"
ax2.yticklabelsvisible = false
ax.yticks = ([0.001,0.002, 0.005, 0.01,0.02, 0.05, 0.1, 0.2, 0.5, 1], ["0.1","0.2", "0.5", "1","2", "5", "10","20", "50","100"])



ii = 14
# Add a rotated label at the chosen point
xx = x



xlims!(ax, 0, 3000)
xlims!(ax2, 0, 3000)
ylims!(ax, 0.001, 1)
ylims!(ax2, 0.001, 1)
ax.xticks = ([30.0,300,600,1000,1500,2000,2550,3000], ["30","300","600","1000","1500","2000","5000","7500"])
ax2.xticks = ([30.0,300,600,1000,1500,2000,2550,3000], ["30","300","600","1000","1500","2000","5000","7500"])

ax3 = Axis(fig[2, 1])
ax4 = Axis(fig[2, 2])
lines!(ax3, x2, data_dry*100, linewidth=3,alpha=1,color=:black)
lines!(ax4, x2, data_wet*100, linewidth=3,alpha=1,color=:black)

rowsize!(fig.layout, 2, Relative(0.25))
xlims!(ax3, 0, 3000)
xlims!(ax4, 0, 3000)
ylims!(ax3, 0.0, 100)
ylims!(ax4, 0.0, 100)
ax3.xticks = ([30.0,300,600,1000,1500,2000,2550,3000], ["30","300","600","1000","1500","2000","5000","7500"])
ax4.xticks = ([30.0,300,600,1000,1500,2000,2550,3000], ["30","300","600","1000","1500","2000","5000","7500"])
ax3.xlabel = "Footprint size (m)"
ax4.xlabel = "Footprint size (m)"
ax.xticklabelsvisible = false
ax2.xticklabelsvisible = false

ax3.yticks = ([0,25,50,75,100], ["0","25","50","75","100"])
ax4.yticks = ([0,25,50,75,100], ["0","25","50","75","100"])
ax4.yticklabelsvisible = false
ax3.ylabel = "Fraction of scenes < 0.1%"
#colsize!(fig.layout, 2, Relative(0.47))
fig
CairoMakie.save("Cloud_WetDrySeason_S2_v3.pdf", fig)

function compAngle(xx,yy,index)
    x1 = xx[index-1]
    x2 = xx[index]
    y1 = yy[index-1]
    y2 = yy[index]
    return atan((y2-y1)/(x2-x1))
end