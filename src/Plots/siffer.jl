using CairoMakie, NCDatasets, GeoMakie
#Reading in TROPOMI data
ncfile =Dataset("/net/fluo/data1/ftp/data/tropomi/gridded/SIF740/TROPOMI-SIF740nm_01-2020--12-2020_0_2deg_1-monthly.nc")
trop_lat = ncfile["lat"][:];
trop_lon = ncfile["lon"][:];
trop_SIF = ncfile["sif"][:];
# Just pick winter months, DJF
trop_SIF = trop_SIF[[1,2,12],:,:];

mean_trop_SIF=zeros(length(trop_lon), length(trop_lat));
Nsamples = zeros(length(trop_lon), length(trop_lat));
for ilon=1:length(trop_lon)
    for ilat=1:length(trop_lat)
        if (trop_lat[ilat]<-60.0)
            for it=1:size(trop_SIF)[1]
                if (!ismissing(trop_SIF[it, ilon, ilat]))
                    mean_trop_SIF[ilon, ilat] += trop_SIF[it, ilon, ilat];
                    Nsamples[ilon, ilat] += 1;
                end
            end
        end
    end
end
for ilon=1:length(trop_lon)
    for ilat=1:length(trop_lat)
        if (Nsamples[ilon,ilat]>=1)
            mean_trop_SIF[ilon,ilat]/=Nsamples[ilon,ilat]
        else
            mean_trop_SIF[ilon,ilat] = NaN
        end
    end
end
#wo = findall(x -> x >= -900, trop_SIF)
#wox = [index[2] for index in wo]
#woy = [index[3] for index in wo]
x = convert(Vector{Float32}, trop_lon)
y = convert(Vector{Float32}, trop_lat)
min_ = minimum(mean_trop_SIF)
max_ = maximum(mean_trop_SIF)
min2 = -0.2 #-0.5*(abs(min_)+abs(max_))#maximum([abs(min), abs(max)]) #0.5*(min+max)
max2 = 0.2 # 0.5*(abs(min_)+abs(max_))#maximum([abs(min), abs(max)]) #0.5*(min+max)
#min2 =-0.5*((min_)+(max_))#maximum([abs(min), abs(max)]) #0.5*(min+max)
#max2 = 0.5*((min_)+(max_))
cmap = :viridis
z = mean_trop_SIF
ind = findall(trop_lat.<=60)
f2 = Figure()
ax22 = GeoAxis(f2[7:9, 3:4],
    title = "TROPOMI SIF (740 nm)",
    aspect = AxisAspect(1),
    dest="+proj=sterea +lat_0=-90")
lines!(ax22, GeoMakie.coastlines()) # plot coastlines from Natural Earth as a reference
xlims!(ax22, -180, 180)
ylims!(ax22, -90, -60)
hm2 = heatmap!(ax22, trop_lon, trop_lat[ind], mean_trop_SIF[:,ind], colormap=cmap,  colorrange=(min2, max2))
cb2 = Colorbar(f2[10,3:4], hm2,  label = L"$\overline{I}_\mathrm{SIF, TOA}\,\,[\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mu\mathrm{m}]$",  labelpadding = 10, ticklabelpad = 5, width = 200, colorrange=(min2, max2), vertical=false)