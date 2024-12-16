# Test:
using Glob, Dates, NCDatasets, CairoMakie, GeoMakie, GeoJSON, Downloads, Interpolations

coasts = GeoJSON.read(read(Downloads.download("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_10m_coastline.geojson"), String))
#counties = GeoJSON.read(read(Downloads.download("https://gist.githubusercontent.com/sdwfrost/d1c73f91dd9d175998ed166eb216994a/raw/e89c35f308cee7e2e5a784e1d3afc5d449e9e4bb/counties.geojson"), String))

# Coordinates for Pasadena, CA
lat_pasadena = 34.1478
lon_pasadena = -118.1445


vulcan = Dataset("/home/cfranken/code/gitHub/CarbonI/data/Vulcan_v3_US_annual_1km_total_mn.nc4")
lat = vulcan["lat"][:]
lon = vulcan["lon"][:]
carbon_emissions = vulcan["carbon_emissions"][:,:,:]
emi = sum(carbon_emissions, dims=3)[:,:,1]
ind = findall(lat .> 33.3 .&& lat .< 34.8 .&& lon .> -118.8 .&& lon .< -117.5)

#lats = minimum(lat):0.01:maximum(lat)
lats = 33.68:0.01:34.57
#lons = minimum(lon):0.01:maximum(lon)
lons = -118.703:0.01:-117.587
carbon_emissions_inter = zeros(length(lats), length(lons)) * NaN
for i in eachindex(lats)
    for j in eachindex(lons)
        indi = findall(abs.(lat .-lats[i]) .< 0.007 .&& abs.(lon .-lons[j]) .< 0.007 )
        if length(indi) > 0
            #@show indi[1], i, j
            carbon_emissions_inter[i,j] = emi[indi[1]]
        end
    end
end




fig = Figure(title="Test")
ga = GeoAxis(fig[1, 1]; dest = "+proj=merc",limits=((-118.703, -117.587), (33.68, 34.57)), yticks=33.6:0.2:34.6,xticks=-118.7:0.2:117.6, ytickformat = "{:.1f}ms",xtickformat = "{:.1f}" ) # dest = "+proj=merc", 
    #ga.xtickformat = x -> @sprintf("%.1f", x)  # Format longitude ticks with one decimal place
    #ga.ytickformat = y -> @sprintf("%.1f", y)  # Format latitude ticks with one decimal place
    #lines!(ga, GeoMakie.coastlines(); strokewidth = 0.7, color=:gold, rasterize = 5)
    xco2 = GeoMakie.heatmap!(ga,  lons,lats, carbon_emissions_inter', colormap = :Reds_9, colorrange=(0, 1e5))
    lines!(ga, coasts; color= :grey, xautolimits=false, yautolimits=false, rasterize=5)
    # Add a star marker for Pasadena
    scatter!(ga,[lon_pasadena], [lat_pasadena], marker = :star5, markersize = 10, color = :blue)
    indi = findall(carbon_emissions_inter' .> 1e6)
    scatter!(ga,lons[[ii[1] for ii in indi]], lats[[ii[2] for ii in indi]], marker = :star6, markersize = carbon_emissions_inter'[indi]/3e5, color = :green)
    #hidexdecorations!(ga, grid=false)
    #hideydecorations!(ga, grid=false)
    #Colorbar(fig[1, 2], xco2, label="WRF-Chem XCO₂ (ppm)", height = Relative(0.85))
    #poly!(ga, counties; color= :grey, xautolimits=false, yautolimits=false, rasterize=5)
    #poly!(ga, GeoMakie.to_multipoly(coasts.geometry); strokewidth = 0.7, color=:gold, rasterize = 5)
    #poly!(ga, GeoMakie.to_multipoly(lakes.geometry); strokewidth = 0.7, color=:blue, rasterize = 5,  xautolimits=false, yautolimits=false)
    fig