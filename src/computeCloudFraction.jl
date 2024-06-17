using GeoArrays, StatsBase, CairoMakie
using Glob
la = GeoArrays.read("/net/fluo/data1/data/Landsat/LandSat_mamiraua_TOA_B7.tif")
files= glob("LandSat_mamiraua_2015*.tif","/net/fluo/data1/data/Landsat/")

function computeCloudFraction(la, block)
    countCloudy    = 0
    countCloudFree = 0
    countTotal = 0
    dim = size(la)
    for i1 in 1:block:dim[1]-block, i2 in 1:block:dim[2]-block
        if ~isnan(la[i1,i2,1])
            countTotal += 1
            if any(la[i1:i1+(block-1),i2:i2+(block-1),1] .> 0.07)
                countCloudy += 1
            else
                countCloudFree += 1
            end
        end 
    end
    return countCloudy, countCloudFree, countTotal
end



cf_all = zeros(100,183,5)
for iFile in eachindex(files)
    @show files[iFile]
    la = GeoArrays.read(files[iFile])
        for i = 1:100, j=1:183
            c, cf, ct = computeCloudFraction(la[:,:,j], i)
            cf_all[i,j,iFile] = cf/ct;
        end
end
for i = 1:100, j=1:183
    c, cf, ct = computeCloudFraction(la[:,:,j], i)
    cf_all[i,j] = cf/ct;
end
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

months = [parse(Float64,split(s,"-")[2]) for s in dd]
djf = findall(months.==12 .|| months.==1 .|| months.==2)
mam = findall(months.==3 .|| months.==4 .|| months.==5)
jja = findall(months.==6 .|| months.==7 .|| months.==8)
son = findall(months.==9 .|| months.==10 .|| months.==11)

cf_median = [percentile(cf_all[i,:],50) for i=1:100] 
cf_20 = [percentile(cf_all[i,:,1],20) for i=1:100]
cf_80 = [percentile(cf_all[i,:,1],80) for i=1:100]

cf_median_djf = [percentile(cf_all[i,djf],50) for i=1:100]
cf_median_mam = [percentile(cf_all[i,mam],50) for i=1:100]
cf_median_jja = [percentile(cf_all[i,jja],50) for i=1:100]
cf_median_son = [percentile(cf_all[i,son],50) for i=1:100]

cf_median_djf = [percentile(cf_all[i,ind[djf],1],50) for i=1:100]
cf_median_mam = [percentile(cf_all[i,ind[mam],1],50) for i=1:100]
cf_median_jja = [percentile(cf_all[i,ind[jja],1],50) for i=1:100]
cf_median_son = [percentile(cf_all[i,ind[son],1],50) for i=1:100]

c = to_color(:gray)
c2 = RGBAf(c.r, c.g, c.b,0.2)
f = Figure(resolution=(400,300))
ax1 = Axis(f[1,1], title="Landsat 8, confidently cloud free statistics", ylabel="Cloud-free pixels (%)", xlabel="Footprint size (m)")

#lines!(ax1, collect(30:30:3000), cf_all[:,10]*100, linewidth=4, alpha=0.7,color="black")
lowery = cf_20*100
uppery = cf_80*100
band!(ax1, collect(30:30:3000), lowery, uppery, transparency=true, fillalpha=0.3, color=c2)

#lines!(ax1, collect(30:30:3000), median(cf_all,dims=2)[:,1]*100,color="black", linewidth=3, alpha=0.7)
lines!(ax1, collect(30:30:3000), cf_median_djf*100, linewidth=3, alpha=0.7, label="DJF")
lines!(ax1, collect(30:30:3000), cf_median_mam*100, linewidth=3, alpha=0.7, label="MAM")
lines!(ax1, collect(30:30:3000), cf_median_jja*100, linewidth=3, alpha=0.7, label="JJA")
lines!(ax1, collect(30:30:3000), cf_median_son*100, linewidth=3, alpha=0.7, label="SON")
lines!(ax1, [60,60], [0,39], linewidth=4, alpha=0.2, color=:black)#, label="Carbon-I Target mode")
lines!(ax1, [300,300], [0,23], linewidth=4, alpha=0.2, color=:black)#, label="Carbon-I Global mode")
lines!(ax1, [1300,1300], [0,3], linewidth=4, alpha=0.2, color=:black)#, label="OCO-2")
#inset_ax2 = add_axis_inset(f[1, 1]; bgcolor=(:white, 0.85),
#halign=:right, valign=:center,
#        width=Relative(0.35), height=Relative(0.3))
#        lines!(inset_ax2, collect(30:30:3000), cf_median_djf*100, linewidth=3, alpha=0.7, label="DJF")
#        lines!(inset_ax2, collect(30:30:3000), cf_median_mam*100, linewidth=3, alpha=0.7, label="MAM")
#        lines!(inset_ax2, collect(30:30:3000), cf_median_jja*100, linewidth=3, alpha=0.7, label="JJA")
#        lines!(inset_ax2, collect(30:30:3000), cf_median_son*100, linewidth=3, alpha=0.7, label="SON")
#        CairoMakie.xlims!(inset_ax2,1200,2090)
#CairoMakie.ylims!(inset_ax2,0,2)

CairoMakie.xlims!(ax1,0,2090)
CairoMakie.ylims!(ax1,0,50)
#xlims!(ax_,2040,2380)
axislegend(ax1, position=:rt)
f

c = to_color(:gray)
cb = to_color(:black)
c2 = RGBAf(c.r, c.g, c.b,0.2)
c3 = RGBAf(cb.r, cb.g, cb.b,0.5)

function add_axis_inset(pos=fig[1, 1]; bgcolor=:snow2,
    halign, valign, width=Relative(0.5),height=Relative(0.35),
    alignmode=Mixed(left=5, right=5))
    inset_box = Axis(pos; width, height, halign, valign, alignmode,
        xticklabelsize=12, yticklabelsize=12, backgroundcolor=bgcolor)
    # bring content upfront
    CairoMakie.translate!(inset_box.scene, 0, 0, 10)
    return inset_box
end