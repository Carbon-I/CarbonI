
using JSON, CairoMakie, MakieThemes, JLD2, Statistics, StatsBase

load("cf_all.jld2")
a = load("cf_all.jld2")
cf_all = a["cf_all"]

j = JSON.parsefile("2015_2019_dates.json")
[
  "2020-01-01",
  "2020-01-17",
  "2020-02-02",
  "2020-02-18",
  "2020-03-21",
  "2020-04-22",
  "2020-05-24",
  "2020-06-09",
  "2020-07-11",
  "2020-07-27",
  "2020-08-12",
  "2020-08-28",
  "2020-09-13",
  "2020-09-29",
  "2020-10-15",
  "2020-10-31",
  "2020-11-16",
  "2020-12-02",
  "2020-12-18",
  "2020-06-25",
  "2020-01-08",
  "2020-01-24",
  "2020-02-09",
  "2020-02-25",
  "2020-03-12",
  "2020-03-28",
  "2020-04-29",
  "2020-05-31",
  "2020-07-02",
  "2020-07-18",
  "2020-08-03",
  "2020-08-19",
  "2020-09-04",
  "2020-09-20",
  "2020-10-06",
  "2020-10-22",
  "2020-11-23",
  "2020-12-09",
  "2020-12-25",
  "2020-04-13"
]

months = [parse(Float64,split(j["DateTime"][s],"-")[2]) for s in keys(j["DateTime"])]
ind    = [parse(Int,s) for s in keys(j["DateTime"])].+1

#months = [parse(Float64,split(s,"-")[2]) for s in dd]
wet = findall(months.==12 .|| months.==1 .|| months.==2 .|| months.==3 .|| months.==4 .|| months.==5 .|| months.==9 .|| months.==10 .|| months.==11)
#mam = findall(months.==3 .|| months.==4 .|| months.==5)
dry = findall(months.==6 .|| months.==7 .|| months.==8)
#son = findall(months.==9 .|| months.==10 .|| months.==11)

cf_median = [percentile(cf_all[i,:,1],50) for i=1:100] 
cf_20 = [percentile(cf_all[i,:,1],20) for i=1:100]
cf_80 = [percentile(cf_all[i,:,1],80) for i=1:100]

cf_median_dry = [percentile(cf_all[i,ind[dry],1][:],50) for i=1:100]
cf_median_wet = [percentile(cf_all[i,ind[wet],1][:],50) for i=1:100]

c = to_color(:gray)
c2 = RGBAf(c.r, c.g, c.b,0.2)

cc = to_color(:black)
c3 = RGBAf(cc.r, cc.g, cc.b,0.7)

Makie.set_theme!(ggthemr(:earth))
f = Figure(resolution=(600,430))
ax1 = Axis(f[1,1], halign=:left, title="Landsat 8, confidently cloud free statistics", ylabel="Cloud-free pixels (%)", xlabel="Footprint size (m)")

#lines!(ax1, collect(30:30:3000), cf_all[:,10]*100, linewidth=4, alpha=0.7,color="black")
lowery = cf_20*100
uppery = cf_80*100
band!(ax1, collect(30:30:3000), lowery, uppery, transparency=true, fillalpha=0.3, color=c2, label="20-80%ile")

#lines!(ax1, collect(30:30:3000), median(cf_all,dims=2)[:,1]*100,color="black", linewidth=3, alpha=0.7)
lines!(ax1, collect(30:30:3000), cf_median_dry*100, linewidth=3, alpha=0.7, label="Dry season")
lines!(ax1, collect(30:30:3000), cf_median_wet*100, linewidth=3, alpha=0.7, label="Wet season")
#lines!(ax1, collect(30:30:3000), cf_median_jja*100, linewidth=3, alpha=0.7, label="JJA")
#lines!(ax1, collect(30:30:3000), cf_median_son*100, linewidth=3, alpha=0.7, label="SON")
#lines!(ax1, [60,60], [0,39], linewidth=4, alpha=0.2, color=c3)#, label="Carbon-I Target mode")
#lines!(ax1, [300,300], [0,23], linewidth=4, alpha=0.2, color=c3)#, label="Carbon-I Global mode")
#lines!(ax1, [1300,1300], [0,3], linewidth=4, alpha=0.2, color=c3)#, label="OCO-2")
#inset_ax2 = add_axis_inset(f[1, 1]; bgcolor=(:white, 0.85),
#halign=:right, valign=:center,
#        width=Relative(0.35), height=Relative(0.3))
#        lines!(inset_ax2, collect(30:30:3000), cf_median_djf*100, linewidth=3, alpha=0.7, label="DJF")
#        lines!(inset_ax2, collect(30:30:3000), cf_median_mam*100, linewidth=3, alpha=0.7, label="MAM")
#        lines!(inset_ax2, collect(30:30:3000), cf_median_jja*100, linewidth=3, alpha=0.7, label="JJA")
#        lines!(inset_ax2, collect(30:30:3000), cf_median_son*100, linewidth=3, alpha=0.7, label="SON")
#        CairoMakie.xlims!(inset_ax2,1200,2090)
#CairoMakie.ylims!(inset_ax2,0,2)

CairoMakie.xlims!(ax1,0,3050)
CairoMakie.ylims!(ax1,0,50)
#xlims!(ax_,2040,2380)
axislegend(ax1, position=:rt)

f
CairoMakie.save("CarbonI_Clouds_Earth.pdf", f)

res = collect(30:30:3000);
res = collect(30:30:3000)
swath_width = res # (*1000/1000) in km (1000 across track pixels)
revisit = 40000 ./(swath_width*15);
# Compute total number of pixels in a 2 degree box:
# 2 degree box is 222.39 km across
nPix = 222.39e3.^2 ./res.^2; ./ swath_width

f = Figure(resolution=(600,430))
ax1 = Axis(f[1,1], title="Theoretical Revisit times", ylabel="Revisit (days)", xlabel="Footprint size (m)", yscale=Makie.pseudolog10, yticks = [1, 7, 30, 60, 100])
ax2 = Axis(f[1,1], xlabel = "x axis", ylabel = "Swath width (km)")
lines!(ax1, res, revisit,  color=:white, linewidth=3, alpha=0.7, label="Revisit")
lines!(ax2, res, swath_width ,  color=:red, linewidth=3, alpha=0.7,label="Swath width", grid = false) 
#lines!(ax1, res, revisit  ,  linewidth=3, alpha=0.7, label="Wet season")
axislegend(ax1, position=:lt)
axislegend(ax2, position=:rb)
CairoMakie.xlims!(ax1,0,3050)
CairoMakie.xlims!(ax2,0,3050)
CairoMakie.ylims!(ax1,1,100)
CairoMakie.ylims!(ax2,0,3000)
ax2.yaxisposition = :right
ax2.yticklabelalign = (:left, :center)
ax2.xticklabelsvisible = false
ax2.xticklabelsvisible = false
ax2.xlabelvisible = false
ax2.ygridvisible = false
#hidedecorations!(ax2)
#hidespines!(ax2)
#hideydecorations!(ax2)
#linkxaxes!(ax1,ax2)
f
CairoMakie.save("CarbonI_TheoreticalRevisit.pdf", f)

f2 = Figure(resolution=(600,430))
ax1 = Axis(f2[1,1], title="Effective Revisit times", ylabel="Revisit (days)", xlabel="Footprint size (m)", yscale=Makie.pseudolog10, yticks = [10, 40, 60, 100, 200, 365])
lines!(ax1, res, revisit ./cf_median_dry , linewidth=3, alpha=0.7, label="Dry season")
lines!(ax1, res, revisit ./cf_median_wet ,  linewidth=3, alpha=0.7, label="Wet season")
axislegend(ax1, position=:lt)
CairoMakie.xlims!(ax1,0,3050)
f2
CairoMakie.save("CarbonI_EffectiveRevisit.pdf", f2)

f2 = Figure(resolution=(600,430))
ax1 = Axis(f2[1,1], title="Effective Revisit times", ylabel="Revisit (days)", xlabel="Footprint size (m)", yscale=Makie.pseudolog10, yticks = [10, 40, 60, 100, 200, 365])
lines!(ax1, res, revisit ./((1 .-cf_median_dry).^(nPix./1000)) , linewidth=3, alpha=0.7, label="Dry season")
lines!(ax1, res, revisit ./((1 .-cf_median_wet).^(nPix./1000)),  linewidth=3, alpha=0.7, label="Wet season")
axislegend(ax1, position=:lt)
CairoMakie.xlims!(ax1,0,3050)
f2
CairoMakie.save("CarbonI_EffectiveRevisit.pdf", f2)

plot(res, cf_median_wet*100 ./revisit)
plot!(res, cf_median_dry*100 ./revisit)
plot!(res, cf_20*100 ./revisit)
plot!(res, cf_80*100 ./revisit)

