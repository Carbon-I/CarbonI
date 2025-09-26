using NCDatasets
using CairoMakie, GeoMakie
using FileIO, Printf

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
AMAZON_LON = (-82.0, -33.0)
AMAZON_LAT = (-13.0,  13.0)

# Colormaps / fixed color ranges
CMAP_CH4 = :vik25
CMAP_CO2 = :vik25
CMAP_CO  = :vik25
CRANGE_CH4 = (1800.0, 2100.0)   # ppb
CRANGE_CO2 = (415.0,  435.0)    # ppm
CRANGE_CO  = (5e-6,    2e-3)    # whatever units your CO field has

FIG_SIZE = (1400, 900)
OUTDIR   = "/home/cfranken/data/frames_amazon"
mkpath(OUTDIR)

# =========================
# LOAD COORDS & PREP INDICES
# =========================
ds  = NCDataset(NCFILE, "r")
lat = collect(ds[VAR_LAT][:])         # may be descending
lon = collect(ds[VAR_LON][:])         # may be 0..360
lon = ((lon .+ 180) .% 360) .- 180    # shift to [-180, 180)

# find indices inside box
ilat = findall(x -> AMAZON_LAT[1] ≤ x ≤ AMAZON_LAT[2], lat)
ilon = findall(x -> AMAZON_LON[1] ≤ x ≤ AMAZON_LON[2], lon)

# sort lat/lon ascending and keep permutation for data
lat_ord = sortperm(lat[ilat])
lon_ord = sortperm(lon[ilon])
ilat = ilat[lat_ord]
ilon = ilon[lon_ord]

lat_sub = lat[ilat]
lon_sub = lon[ilon]

# =========================
# FRAME RENDERER
# =========================
function save_frame(d::Int, h::Int; frame::Int)
    # Slice numeric fields (dims: lat, lon, day, hour)
    @show  d, h
    CH4 = ds[VAR_CH4][ilon, ilat, d, h]
    CO2 = ds[VAR_CO2][ilon, ilat, d, h]
    CO  = ds[VAR_CO][ ilon, ilat, d, h]

    fig = Figure(size = FIG_SIZE)

    # Layout: big CH4 left, CO2/CO stacked right
    ga  = GeoAxis(fig[1:2, 1];
        limits = (AMAZON_LON..., AMAZON_LAT...),
        title = "CH₄ over Amazon — day=$(d), slot=$(h)",
        aspect = DataAspect()
    )
    ga2 = GeoAxis(fig[1, 2];
        limits = (AMAZON_LON..., AMAZON_LAT...),
        title = "CO₂ (zoom)",
        aspect = DataAspect()
    )
    ga3 = GeoAxis(fig[2, 2];
        limits = (AMAZON_LON..., AMAZON_LAT...),
        title = "CO (zoom)",
        aspect = DataAspect()
    )

    #colsize!(fig.layout, 1, Relative(0.75))
    #colsize!(fig.layout, 2, Relative(0.25))

    heatmap!(ga,  lon_sub, lat_sub, CH4; interpolate=false,
             colormap=CMAP_CH4, colorrange=CRANGE_CH4)
    heatmap!(ga2, lon_sub, lat_sub, CO2; interpolate=false,
             colormap=CMAP_CO2, colorrange=CRANGE_CO2)
    heatmap!(ga3, lon_sub, lat_sub, CO;  interpolate=false,
             colormap=CMAP_CO,  colorrange=CRANGE_CO)

    # optional coastlines overlay
    for g in (ga, ga2, ga3)
        lines!(g, GeoMakie.coastlines(), linewidth=1.0, color=:black)
    end

    # save with 1D frame counter
    fpath = joinpath(OUTDIR, @sprintf("frame_%03d.png", frame))
    save(fpath, fig)
end

# =========================
# RENDER ALL FRAMES
# =========================
frame = 0
for d in DAYS, h in HOURS
    frame += 1
    @info "Rendering frame=$frame (d=$d, slot=$h)"
    save_frame(d, h; frame=frame)
end

#close(ds)

# After rendering, make a movie with:
# ffmpeg -framerate 4 -i 'frames_amazon/frame_%03d.png' -pix_fmt yuv420p amazon.mp4