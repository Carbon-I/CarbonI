
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
AMAZON_LON = (65.0, 125.0)
AMAZON_LAT = (5.0,  40.0)

# Colormaps / fixed color ranges
CMAP_CH4 = :vik25
CMAP_CO2 = :vik25
CMAP_CO  = :vik25
CRANGE_CH4 = (1800.0, 2150.0)   # ppb
CRANGE_CO2 = (415.0,  440.0)    # ppm
CRANGE_CO  = (0.25,    3.)    # whatever units your CO field has



# Zoom sub-area (inside Amazon bounds)
ZOOM_LON = AMAZON_LON
ZOOM_LAT = AMAZON_LAT


OUTDIR   = "/home/cfranken/data/frames_amazon"
mkpath(OUTDIR)

# =========================
# LOAD DATA & PREP AXES
# =========================
ds  = Dataset(NCFILE, "r")
lat = collect(ds[VAR_LAT][:])    # may be descending (90..-90)
lon = collect(ds[VAR_LON][:])    # may be 0..360
lon = ((lon .+ 180) .% 360) .- 180
time = collect(ds[VAR_TIME][:])  # not strictly used here




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
CO  = ds["tcco"][ ilon, ilat, T1, T2].*1e3

# Figure & tight layout (left wide, right slim)

# Colormaps
CMAP_CH4 = cgrad(:vik25, 11; categorical=true)
CMAP_CO2 = cgrad(:viridis, 11; categorical=true)
CMAP_CO  = cgrad(:Oranges, 11; categorical=true)


using CairoMakie, GeoMakie

function save_frame(d::Int, h::Int; frame::Int)
    # Slice numeric fields (dims: lat, lon, day, hour)
    @show  d, h
    CH4 = ds[VAR_CH4][ilon, ilat, d, h]
    CO2 = ds[VAR_CO2][ilon, ilat, d, h]
    CO  = ds[VAR_CO][ ilon, ilat, d, h].*1000

    fig = plotAmazon(CH4, CO2, CO)

    

    # save with 1D frame counter
    fpath = joinpath(OUTDIR, @sprintf("asia_frame_%03d.png", frame))
    save(fpath, fig; dpi=300)
end


function plotAmazon(CH4, CO2, CO)
  @show maximum(CH4)
  fig = Figure(resolution = (1000, 400))

  # Nest a layout for the right column and split it 50/50
  rightcol = GridLayout(fig[1:2, 2])
  rowsize!(rightcol, 1, Relative(0.5))
  #rowsize!(rightcol, 2, Relative(0.5))
  #rowsize!(rightcol, 2, Relative(0.5))
  rowgap!(rightcol, 0)
  colgap!(rightcol, 0)

  colsize!(fig.layout, 1, Relative(0.65))
  colsize!(fig.layout, 2, Relative(0.35))
  rowgap!(fig.layout, 0)
  colgap!(fig.layout, 0)

  # Common GeoAxis settings: rectangular lon/lat, no decorations, tight box
  common_geo = (
      source = "+proj=latlong +datum=WGS84",
      dest   = "+proj=eqc",
      limits = (AMAZON_LON..., AMAZON_LAT...),
      aspect = DataAspect(),
  )

  # Left axis spans both rows
  ga  = GeoAxis(fig[1:2, 1]; common_geo...)
  hidedecorations!(ga); hidespines!(ga)

  # Right-top / right-bottom axes live in the nested layout
  ga2 = GeoAxis(rightcol[1, 1]; common_geo...)
  ga3 = GeoAxis(rightcol[2, 1]; common_geo...)
  for g in (ga2, ga3)
      hidedecorations!(g); hidespines!(g)
  end

  # ---- Data layers ----
  hm1 = heatmap!(ga,  lon_sub, lat_sub, CH4; interpolate=false, colormap=CMAP_CH4, colorrange=CRANGE_CH4)
  hm2 = heatmap!(ga2, lon_sub, lat_sub, CO2; interpolate=false, colormap=CMAP_CO2, colorrange=CRANGE_CO2)
  hm3 = heatmap!(ga3, lon_sub, lat_sub, CO;  interpolate=false, colormap=CMAP_CO,  colorrange=CRANGE_CO)

  tightlimits!(ga); tightlimits!(ga2); tightlimits!(ga3)
  

  # Coastlines (don’t change limits)
  for g in (ga, ga2, ga3)
      lines!(g, GeoMakie.coastlines(); color=:black, linewidth=1.8,
            xautolimits=false, yautolimits=false)
  end

  

  # ---- Titles INSIDE each map so they don't steal layout height ----
  # (avoid Label(fig[*,*,Top()]) which creates extra row height)
  text!(ga,  AMAZON_LON[2] - 15, AMAZON_LAT[2] - 1, text="CH₄",  align=(:left,:top),   color=:black, font = :bold, fontsize=18)
  text!(ga2, AMAZON_LON[2] - 16, AMAZON_LAT[2] - 1, text="CO₂",  align=(:left,:top),   color=:black, font = :bold, fontsize=16)
  text!(ga3, AMAZON_LON[2] - 16, AMAZON_LAT[2] - 1, text="CO",   align=(:left,:top),   color=:black, font = :bold, fontsize=16)

  colsize!(fig.layout, 1, Relative(0.7))
  colsize!(fig.layout, 2, Relative(0.3))
  rowgap!(fig.layout, 0); colgap!(fig.layout, 0)
  fig
end


frame = 0
for d in DAYS, h in HOURS
    frame += 1
    @info "Rendering frame=$frame (d=$d, slot=$h)"
    save_frame(d, h; frame=frame)
end