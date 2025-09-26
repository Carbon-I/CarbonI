
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
AMAZON_LON = (-82.0, -33.0)
AMAZON_LAT = (-15.0,  13.0)

# Colormaps / fixed color ranges
CMAP_CH4 = :vik25
CMAP_CO2 = :vik25
CMAP_CO  = :vik25
CRANGE_CH4 = (1800.0, 2100.0)   # ppb
CRANGE_CO2 = (415.0,  435.0)    # ppm
CRANGE_CO  = (5e-6,    2e-3)    # whatever units your CO field has

# Amazon basin bounds (degrees)
AMAZON_LON = (-82.0, -33.0)
AMAZON_LAT = (-13.0,  13.0)

# Zoom sub-area (inside Amazon bounds)
ZOOM_LON = (-82.0, -33.0)
ZOOM_LAT = (-13.0,  13.0)





# =========================
# LOAD DATA & PREP AXES
# =========================
ds  = Dataset(NCFILE, "r")
lat = collect(ds[VAR_LAT][:])    # may be descending (90..-90)
lon = collect(ds[VAR_LON][:])    # may be 0..360
lon = ((lon .+ 180) .% 360) .- 180
time = collect(ds[VAR_TIME][:])  # not strictly used here



# --- choose your bounding box ---
AMAZON_LON = (-82.0 , -33.0 )
AMAZON_LAT = (-14.0,  13.0)

# --- get the index ranges for Amazon only ---
# find indices where coords fall inside the box
ilat = findall(x -> AMAZON_LAT[1] ≤ x ≤ AMAZON_LAT[2], lat)
ilon = findall(x -> AMAZON_LON[1] ≤ x ≤ AMAZON_LON[2], lon)

# crop arrays
lat_sub = lat[ilat]
lon_sub = lon[ilon]

# --- now slice only the subwindow for one timestep ---
T1, T2 = 10, 6
CH4 = ds["tcch4"][ilon, ilat, T1, T2]
CO2 = ds["tcco2"][ilon, ilat, T1, T2]
CO  = ds["tcco"][ ilon, ilat, T1, T2]

# Figure & tight layout (left wide, right slim)

# Colormaps
CMAP_CH4 = cgrad(:vik25, 11; categorical=true)
CMAP_CO2 = cgrad(:vik25, 11; categorical=true)
CMAP_CO  = cgrad(:Oranges, 11; categorical=true)

# Fixed value ranges (set to `nothing` to auto from data windows below)
CRANGE_CH4 = (1880.0, 2080.0)   # ppb
CRANGE_CO2 = (420.0,  430.0)    # ppm
CRANGE_CO  = ( 0.0005,  0.002)    # ppb


fig = Figure(resolution = (1400, 500))


common_geo_kwargs = (
    source = "+proj=latlong +datum=WGS84",  # data CRS (lon/lat)
    dest   = "+proj=eqc",                   # Plate Carrée
    limits = (AMAZON_LON..., AMAZON_LAT...),
    xticksvisible = false, yticksvisible = false,
    xticklabelsvisible = false, yticklabelsvisible = false,
    xlabelvisible = false, ylabelvisible = false,
    leftspinevisible = false, rightspinevisible = false,
    topspinevisible = false, bottomspinevisible = false,
    aspect = DataAspect(),
)

# Axes
# Axes WITHOUT titles (we'll place titles with Label to avoid extra axis padding)
ga  = GeoAxis(fig[1:2, 1]; common_geo_kwargs...)
ga2 = GeoAxis(fig[1,   2]; common_geo_kwargs...)
ga3 = GeoAxis(fig[2,   2]; common_geo_kwargs...)

# Make the axis box equal to just the plotting area
ga.alignmode  = Inside()
ga2.alignmode = Inside()
ga3.alignmode = Inside()

# Hide all decorations/spines so no extra space is reserved
for g in (ga, ga2, ga3)
    hidedecorations!(g)
    hidespines!(g)
end




# Heatmaps (use your cmaps/clims)
hm1 = heatmap!(ga,  lon_sub, lat_sub, CH4; interpolate=false, colormap=CMAP_CH4, colorrange=CRANGE_CH4)
hm2 = heatmap!(ga2, lon_sub, lat_sub, CO2; interpolate=false, colormap=CMAP_CO2, colorrange=CRANGE_CO2)
hm3 = heatmap!(ga3, lon_sub, lat_sub, CO;  interpolate=false, colormap=CMAP_CO,  colorrange=CRANGE_CO)

tightlimits!(ga); tightlimits!(ga2); tightlimits!(ga3)

# Coastlines overlaid (don’t let them change limits)
for g in (ga, ga2, ga3)
    lines!(g, GeoMakie.coastlines(); color=:black, linewidth=1.8,
           xautolimits=false, yautolimits=false)
end


# Colorbars flush-right of each axis (no gaps)
cb1 = Colorbar(fig[1:2,1], hm1; label="CH₄ (ppb)", width=12, tellheight=false, tellwidth=false)
#cb2 = Colorbar(fig[1,2], hm2; label="CO₂ (ppm)", width=12, tellheight=false, tellwidth=false)
#cb3 = Colorbar(fig[2,2], hm3; label="CO (ppb)",  width=12, tellheight=false, tellwidth=false)

for cb in (cb1,)
  cb.alignmode = Inside()              # position inside the plot area
cb.halign    = :0.9                # right edge
cb.valign    = :0.65                  # top edge
cb.height    = Relative(0.30)        # 30% of the map height (adjust as you like)
cb.width     = 12                    # bar thickness in px
end

for cb in (cb2, cb3)
  cb.alignmode = Inside()              # position inside the plot area
cb.halign    = :0.9                # right edge
cb.valign    = :0.6                  # top edge
cb.height    = Relative(0.30)        # 30% of the map height (adjust as you like)
cb.width     = 12                    # bar thickness in px
end

colsize!(fig.layout, 1, Relative(0.7))
colsize!(fig.layout, 2, Relative(0.3))
rowgap!(fig.layout, 0); colgap!(fig.layout, 0)

fig