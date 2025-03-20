using DelimitedFiles, CairoMakie, Statistics, LaTeXStrings
using CSV, DataFrames, Glob, Statistics,  StatsBase
using Makie, CairoMakie

function extract_month(date_string::String)
    # Extracting the month part from the string
    month_str = date_string[5:6]
    # Converting the month string to integer
    parse(Int, month_str)
end

# Function to check if the month is in the range July to September
is_in_dry(row) = row.Month in 7:9
is_in_wet(row) = row.Month in [1,2,3,4,5,12]

x = [30,90,150, 200, 300, 400, 600, 800, 1000, 1250, 1500, 1750, 2000,3000, 5000, 7500]#, 3000, 5000, 7500, 10000];
datasets = []
for res in (30,90,150, 200, 300, 400, 600, 800, 1000, 1250, 1500, 1750, 2000,3000, 5000, 7500)#, 3000, 5000, 7500, 10000)
#for res in (200,  2000)#, 3000, 5000, 7500, 10000)
    files = Glob.glob("AmazonBox2018-01-01_2022-01-01_*"* string(res)*".0.csv","/home/cfranken/shared/GE_BOX_final/")
    df = CSV.read(files, DataFrame; source = "source")
    df[!,:Month] = [extract_month(date_string) for date_string in df.var"system:index"]

    push!(datasets,df)
end

med_dry = [quantile(filter(is_in_dry,d).cloud_free_fraction02,[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]) for d in datasets]
med_wet = [quantile(filter(is_in_wet,d).cloud_free_fraction02,[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]) for d in datasets]
med_all = [quantile(d.cloud_free_fraction02,[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]) for d in datasets]

cf_mean_dry = [mean(filter(is_in_dry,d).cloud_free_fraction02) for d in datasets]
cf_mean_wet = [mean(filter(is_in_wet,d).cloud_free_fraction02) for d in datasets]

percentiles_wet = mapreduce(permutedims, vcat, med_wet)
percentiles_dry = mapreduce(permutedims, vcat, med_dry)
percentiles_all = mapreduce(permutedims, vcat, med_all )
CloudStats = [x percentiles_all[:,5]];
writedlm("data/CloudStatsCentralAmazonia.dat", CloudStats, delim='\t')

include(joinpath(@__DIR__, "Plots", "CI_colors.jl"))

cs = readdlm("data/tropics_month_res_stats.csv", ',',skipstart=1)

area    = 12^2
gp_area = 0.4^2

npix_cloudFree = area/gp_area

resi = unique(cs[:,2])
months = unique(cs[:,1])

med_cf = zeros(length(resi))
for (ir, res) in enumerate(resi)
    i = findall(x->x==res, cs[:,2])
    med_cf[ir] = mean(cs[i,4])/2
end

# Test wih old data:
med_cf = (percentiles_wet[:,5] + percentiles_dry[:,5])/2
resi = x

sigma_ss = 1.5:0.1:4

sigma_aggregate = sigma_ss' ./ sqrt.(npix_cloudFree .* med_cf)

f = Figure(resolution=(400,500), title="Error in 12km x 12km area (%)", fontsize=14)
ax = Axis(f[1,1], xlabel="GSD (m) (only affecting cloud fraction)", ylabel=L"\text{Median f_{Cloud Free} (%)}", title="Cloud statistics in the central Amazon (worst case)") 
co_ = CairoMakie.lines!(ax, resi, med_cf*100, linewidth = 3, color=CarbonI_colors[1])
CairoMakie.xlims!(resi[1],1000)
CairoMakie.ylims!(0,17)
ax1 = Axis(f[2,1], xlabel=L"\text{GSD (m)}", ylabel=L"\text{\sigma_s (%)}", )

co = CairoMakie.contourf!(ax1, resi, sigma_ss, sigma_aggregate * 1900/100,levels=2:25)#, levels=[0.2, 0.3,0.4, 0.5, 0.6,0.7],  labels=true, colorrange=(0.01,0.5), colormap = (:viridis, 0.5), extendhigh = (:orange,0.4), extendlow = (:gray,0.4)); 
CairoMakie.contour!(ax1, resi, sigma_ss, sigma_aggregate * 1900/100, levels=[8,12],  labels=true,  labelsize = 14,labelfont = :bold, labelcolor = :black, color=:black); 
#S = findall(scales .< 25)
#CairoMakie.contourf!(ax2, scales[iS ], cloud_fracs,1.0./sqrt.(eff_pix_area[iS,:]./pixArea)*n2o_error(0.15)*100, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), colormap = (:viridis, 0.5), extendhigh = (:orange,0.4), extendlow = (:gray,0.4)); 
#CairoMakie.contour!(ax2, scales[iS ], cloud_fracs,1.0./sqrt.(eff_pix_area[iS,:]./pixArea)*n2o_error(0.15)*100, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), labelsize = 14,labelfont = :bold); 

CairoMakie.Colorbar(f[2,2],co,   label=L"\text{CH_4 1\sigma error (ppb) in 12km x 12km area}")
CairoMakie.ylims!(sigma_ss[1],3.99)
CairoMakie.xlims!(resi[1],1000)
x1 = 400
y1 = 3.7
#scatter!(ax1,x1, y1)
text!(x1, y1, text = "Requirement", align = (:left,:bottom),font = :bold,color = :white, rotation=-0.8)

x2 = 300
y2 = 2.7
scatter!(ax1,x2, y2, color=:orange)
text!(x2, y2, text = "CBE", align = (:left,:top),font = :bold,color = :white)
#hideydecorations!(ax2, grid=false)
hidexdecorations!(ax, grid=false)
#hideydecorations!(ax2, grid=false)

f
save("plots/AggregationErrors.pdf", f)
save("plots/AggregationErrors.png", f, dpi=400)



    