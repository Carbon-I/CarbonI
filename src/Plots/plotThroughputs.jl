using DelimitedFiles, CairoMakie, Statistics, LaTeXStrings
using CSV, DataFrames, Glob, Statistics,  StatsBase
using Makie, CairoMakie

include("src/Plots/CI_colors.jl")

T = readdlm("data/throughput_comparison.txt")
#Telescope (3 reflections)	Bandpass Filter	Dyson Lens	Grating	QE	System Throughput w/o QE	Total Throughput
oo = readdlm("data/optics.txt", skipstart=1)

d =0.02
f = Figure(resolution=(450,450), title="Throughput", fontsize=14)
ax0 = Axis(f[1,1], xlabel="Wavelength (nm)", ylabel="Per element", title="Optical Efficiences vs Wavelength")
ax = Axis(f[2,1], xlabel="Wavelength (nm)", ylabel="Total Efficiency", yminorticks = IntervalsBetween(5),yminorgridvisible = true) 
em = CairoMakie.lines!(ax, T[:,1],T[:,2], linewidth = 3, color=CarbonI_colors[10], label="EMIT")
text!(T[10,1],T[10,2]+d, text = "EMIT", align = (:left,:bottom),font = :bold,color = CarbonI_colors[10])
cm = CairoMakie.lines!(ax, T[:,1],T[:,3], linewidth = 3, color=CarbonI_colors[4], label="CarbonMapper")
text!(T[10,1],T[10,3]+d, text = "CarbonMapper", align = (:left,:bottom),font = :bold,color = CarbonI_colors[4])
ci = CairoMakie.lines!(ax, T[:,1],T[:,4], linewidth = 3, color=CarbonI_colors[1], label="Carbon-I")
text!(T[10,1],T[10,4]+d, text = "Carbon-I", align = (:left,:bottom),font = :bold,color = CarbonI_colors[1])

CairoMakie.lines!(ax0, oo[:,1]*1e3,oo[:,2], linewidth = 2, color=CarbonI_colors[2], label="Telescope")
CairoMakie.lines!(ax0, oo[:,1]*1e3,oo[:,3], linewidth = 2, color=CarbonI_colors[5], label="Bandpass")
CairoMakie.lines!(ax0, oo[:,1]*1e3,oo[:,4], linewidth = 2, color=CarbonI_colors[11], label="Dyson")
CairoMakie.lines!(ax0, oo[:,1]*1e3,oo[:,5], linewidth = 2, color=CarbonI_colors[6], label="Grating")
axislegend(ax0, position=:rb, labelsize =10,orientation = :horizontal)
CairoMakie.xlims!(ax,T[1,1],T[end,1])
CairoMakie.xlims!(ax0,T[1,1],T[end,1])
CairoMakie.ylims!(ax,0,1)
CairoMakie.ylims!(ax0,0.8,1)
hidexdecorations!(ax0, grid=false)
f
save("plots/CI-throughputs.pdf", f)
