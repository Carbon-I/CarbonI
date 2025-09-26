using NCDatasets
using CairoMakie, GeoMakie
using Colors, ColorSchemes
using Statistics
using Printf, FileIO

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



# Choose ONE time index to plot
T1 = 10 #(Day, 27 indices)
T2 = 6 #(3hrs step, 8 indices)

# Amazon basin bounds (degrees)
AMAZON_LON = (-82.0, -33.0)
AMAZON_LAT = (-13.0,  13.0)

# Zoom sub-area (inside Amazon bounds)
ZOOM_LON = (-82.0, -33.0)
ZOOM_LAT = (-13.0,  13.0)

# Colormaps
CMAP_CH4 = :vik25
CMAP_CO2 = :vik25
CMAP_CO  = :Oranges

# Fixed value ranges (set to `nothing` to auto from data windows below)
CRANGE_CH4 = (1800.0, 2100.0)   # ppb
CRANGE_CO2 = (415.0,  435.0)    # ppm
CRANGE_CO  = ( 0.00005,  0.002)    # ppb

FIG_SIZE = (1250, 700)  # wide layout



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



# Build index windows
#i_amz, j_amz   = window_indices(lat_asc, lon_asc, AMAZON_LAT[1], AMAZON_LAT[2], AMAZON_LON[1], AMAZON_LON[2])
#i_zoom, j_zoom = window_indices(lat_asc, lon_asc, ZOOM_LAT[1], ZOOM_LAT[2], ZOOM_LON[1], ZOOM_LON[2])

# Extract single-timestep fields
#CH4 = ds[VAR_CH4][:,:,T1,T2]
#CO2 = ds[VAR_CO2][:,:,T1,T2]
#CO  = ds[VAR_CO][:,:,T1,T2]

# =========================
# FIGURE / PLOTS
# =========================
fig = Figure(size = FIG_SIZE)

# Layout: big left (row 1..2, col 1), two stacked right (rows 1 & 2, col 2)
ga  = GeoAxis(fig[1:2, 1];
    limits = (AMAZON_LON..., AMAZON_LAT...),
    title = "CH₄ over Amazon",
    aspect = DataAspect()
)
ga2 = GeoAxis(fig[1, 2];
    limits = (ZOOM_LON..., ZOOM_LAT...),
    title = "CO₂ (zoom)",
    aspect = DataAspect()
)
ga3 = GeoAxis(fig[2, 2];
    limits = (ZOOM_LON..., ZOOM_LAT...),
    title = "CO  (zoom)",
    aspect = DataAspect()
)

# Make columns 58% / 42%
colsize!(fig.layout, 1, Relative(0.4))
colsize!(fig.layout, 2, Relative(0.6))

# Optional coastlines
#GeoMakie.coastlines(ga,  color=(:white, 0.3), linewidth=0.5)
#GeoMakie.coastlines(ga2, color=(:white, 0.3), linewidth=0.5)
#GeoMakie.coastlines(ga3, color=(:white, 0.3), linewidth=0.5)

# Heatmaps with fixed color ranges (numeric → colorbars work)
hm_ch4 = heatmap!(ga,  lon_sub,  lat_sub,  CH4; interpolate=false,
                  colormap=CMAP_CH4, colorrange=CRANGE_CH4)
hm_co2 = heatmap!(ga2, lon_sub, lat_sub, CO2;   interpolate=false,
                  colormap=CMAP_CO2, colorrange=CRANGE_CO2)
hm_co  = heatmap!(ga3, lon_sub, lat_sub, CO;    interpolate=false,
                  colormap=CMAP_CO,  colorrange=CRANGE_CO)

for gg in (ga, ga2, ga3)
    lines!(gg, GeoMakie.coastlines(), linewidth=2.5, color=:black)
end
# Colorbars (inherit colormap from plots; keep labels here)

Colorbar(ga, hm_ch4; label="CH₄ (ppb)", width=12, tellheight=true, tellwidth=false)

cb1 = Colorbar(fig[3, 1], hm_ch4; label = "CH₄ (ppb)")
cb2 = Colorbar(fig[1, 3], hm_co2; label = "CO₂ (ppm)")
cb3 = Colorbar(fig[2, 3], hm_co;  label = "CO (molec/cm2)")

# Tidy colorbar placements (optional: align and shrink)
#rowsize!(fig.layout, 2, Auto(0.2))  # small row for bottom cb
#colsize!(fig.layout, 2, Relative(0.25))  # slim col for right cbs

display(fig)

close(ds)