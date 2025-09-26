
using NCDatasets
using CairoMakie, GeoMakie
using FileIO, Printf
using GridLayoutBase: inset!
# =========================
# USER CONFIG
# =========================
NCFILE   = "/home/cfranken/data/data_sfc_all.nc"
VAR_LAT  = "latitude"
VAR_LON  = "longitude"
VAR_TIME = "valid_time"
VAR_CH4  = "tcch4"
VAR_CO2  = "tcco2"
VAR_CO   = "tcco"

# time indices (27 days, 8 three-hour slots)
DAYS  = 1:27
HOURS = 1:8

# Amazon bounds
AMAZON_LON = (-179, 179)
AMAZON_LAT = (-89.0,  89.0)

# Colormaps / fixed color ranges
CMAP_CH4 = :vik25
CMAP_CO2 = :vik25
CMAP_CO  = :vik25
CRANGE_CH4 = (1850.0, 2150.0)   # ppb
CRANGE_CO2 = (415.0,  435.0)    # ppm
CRANGE_CO  = (0.25,    2.5)    # whatever units your CO field has



# Zoom sub-area (inside Amazon bounds)
ZOOM_LON = AMAZON_LON
ZOOM_LAT = AMAZON_LAT


OUTDIR   = "/home/cfranken/data/frames_amazon"
mkpath(OUTDIR)

# =========================
# LOAD DATA & PREP AXES
# =========================
ds  = Dataset(NCFILE, "r")
lat = Float64.(collect(ds[VAR_LAT][:]))
lon = Float64.(collect(ds[VAR_LON][:]))
lon = ((lon .+ 180) .% 360) .- 180
time = collect(ds[VAR_TIME][:])  # not strictly used here




# --- get the index ranges for Amazon only ---
# find indices where coords fall inside the box
ilat = findall(x -> AMAZON_LAT[1] ≤ x ≤ AMAZON_LAT[2], lat)
ilon = findall(x -> AMAZON_LON[1] ≤ x ≤ AMAZON_LON[2], lon)
step = 2
#ilat = ilat[1:step:end]
#ilon = ilon[1:step:end]
# crop arrays
lat_sub = lat[ilat]
lon_sub = lon[ilon]
p_lon = sortperm(lon_sub)              # sort longitudes ascending
p_lat = sortperm(lat_sub)              # sort latitudes ascending

lon_sorted = lon_sub[p_lon]
lat_sorted = lat_sub[p_lat]

# Figure & tight layout (left wide, right slim)

# Colormaps
CMAP_CH4 = cgrad(:vik25, 21; categorical=true)
CMAP_CO2 = cgrad(:viridis, 21; categorical=true)
CMAP_CO  = cgrad(:vik25, 21; categorical=true)


using CairoMakie, GeoMakie

function save_frame(d::Int, h::Int; frame::Int)
    # Slice numeric fields (dims: lat, lon, day, hour)
    @show  d, h
  
    # CH4 = ds[VAR_CH4][ilon, ilat, d, h]
    #CO2 = ds[VAR_CO2][ilon, ilat, d, h]
    CO_ = ds[VAR_CO][ilon, ilat, d, h]
    CO  = CO_[p_lon, p_lat].*1000


    fig = plotAmazon(CO)

    

    # save with 1D frame counter
    fpath = joinpath(OUTDIR, @sprintf("Global_frame_%03d.png", frame))
    save(fpath, fig; dpi=400)
end


function plotAmazon(CO)
  @show maximum(CO)
  fig = Figure(resolution = (850, 550), backgroundcolor = :black)

  # Left axis spans both rows
  ga  = GeoAxis(fig[1,1], backgroundcolor = :black)
  hidedecorations!(ga); hidespines!(ga)
  hm3 = heatmap!(ga, lon_sorted, lat_sorted, CO;  interpolate=false, colormap=CMAP_CO,  colorrange=CRANGE_CO)



  # Coastlines (don’t change limits)
  for g in (ga,)
      lines!(g, GeoMakie.coastlines(); color=:black, linewidth=1.8,
            xautolimits=false, yautolimits=false)
  end

  fig
end


frame = 0
for d in DAYS, h in HOURS
    frame += 1
    @info "Rendering frame=$frame (d=$d, slot=$h)"
    save_frame(d, h; frame=frame)
end