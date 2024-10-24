# Test:
using Glob, Dates, NCDatasets, CairoMakie, GeoMakie, GeoJSON, Downloads

coasts = GeoJSON.read(read(Downloads.download("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_10m_coastline.geojson"), String))
#counties = GeoJSON.read(read(Downloads.download("https://gist.githubusercontent.com/sdwfrost/d1c73f91dd9d175998ed166eb216994a/raw/e89c35f308cee7e2e5a784e1d3afc5d449e9e4bb/counties.geojson"), String))

# Coordinates for Pasadena, CA
lat_pasadena = 34.1478
lon_pasadena = -118.1445

function extractXCO2(fileNames, outname, time)
    #lats = minimum(lat):0.01:maximum(lat)
    lats = 33.68:0.01:34.57
    #lons = minimum(lon):0.01:maximum(lon)
    lons = -118.703:0.01:-117.587
    xco2_interp = zeros(length(lats), length(lons)) * NaN
    
    for fileName in fileNames
        if isfile(fileName)
            f = Dataset(fileName)
            @show fileName
            co2 = f["CO2_ANT"][:,:,:,1]

            lat = f["XLAT"][:,:,1]  
            lon = f["XLONG"][:,:,1]

            landmask = f["LANDMASK"][:,:,1]
            albedo = f["ALBEDO"][:,:,1]

            levels = f["C3F"][:,1]
            dp = abs.(diff(levels))
            dp = dp/sum(dp)

            xco2 = zeros(size(co2,1), size(co2,2));
            for i in 1:size(co2,1)
                for j in 1:size(co2,2)
                    xco2[i,j] = co2[i,j,:]' *dp;
                end
            end

            

            for i in eachindex(lats)
                for j in eachindex(lons)
                    indi = findall(abs.(lat .-lats[i]) .< 0.007 .&& abs.(lon .-lons[j]) .< 0.007 )
                    if length(indi) > 0
                        #@show indi[1], i, j
                        xco2_interp[i,j] = xco2[indi[1]]
                    end
                end
            end
            #Plots.heatmap!(lons, lats, xco2_interp, clims=(380, 387))
            close(f)
        end
    end
    fig = Figure()
    #lats = 33.68:0.01:34.57
    #lons = minimum(lon):0.01:maximum(lon)
    #lons = -118.703:0.01:-117.587

    ga = GeoAxis(fig[1, 1]; dest = "+proj=merc",limits=((-118.703, -117.587), (33.68, 34.57)), yticks=33.6:0.2:34.6,xticks=-118.7:0.2:117.6, ytickformat = "{:.1f}ms",xtickformat = "{:.1f}" ) # dest = "+proj=merc", 
    #ga.xtickformat = x -> @sprintf("%.1f", x)  # Format longitude ticks with one decimal place
    #ga.ytickformat = y -> @sprintf("%.1f", y)  # Format latitude ticks with one decimal place
    #lines!(ga, GeoMakie.coastlines(); strokewidth = 0.7, color=:gold, rasterize = 5)
    xco2 = GeoMakie.heatmap!(ga,  lons,lats, xco2_interp', colormap = :Oranges, colorrange=(380, 387))
    lines!(ga, coasts; color= :grey, xautolimits=false, yautolimits=false, rasterize=5)
    # Add a star marker for Pasadena
    scatter!(ga,[lon_pasadena], [lat_pasadena], markershape = :star5, markersize = 10, color = :red)
    
    #hidexdecorations!(ga, grid=false)
    #hideydecorations!(ga, grid=false)
    Colorbar(fig[1, 2], xco2, label="WRF-Chem XCO₂ (ppm)", height = Relative(0.85))
    #poly!(ga, counties; color= :grey, xautolimits=false, yautolimits=false, rasterize=5)
    #poly!(ga, GeoMakie.to_multipoly(coasts.geometry); strokewidth = 0.7, color=:gold, rasterize = 5)
    #poly!(ga, GeoMakie.to_multipoly(lakes.geometry); strokewidth = 0.7, color=:blue, rasterize = 5,  xautolimits=false, yautolimits=false)
    fig
    #Plots.title!("XCO2 at $time")
    #savefig(outname)
    save(outname, fig)
end
#=
matching_files = glob("wrfout_d01_2018-01*0401", "/home/cfranken/data/WRF-CO2/")
for file in matching_files
    outname = replace(file, "WRF-CO2/" => "WRF-CO2/plots/")
    outname = outname * ".png"
    @show file, outname
    extractXCO2(file, outname)
end


filenames = [
    "wrfout_d01_2018-01-11_12_00_00_0402", "wrfout_d01_2018-01-11_12_00_00_0427",
    "wrfout_d01_2018-01-13_04_00_00_0402", "wrfout_d01_2018-01-13_04_00_00_0427",
    "wrfout_d01_2018-01-16_13_00_00_0377", "wrfout_d01_2018-01-16_13_00_00_0402"
]
=#
# List of tile digits you are interested in
tile_digits_to_include = ["0402", "0427", "0377", "0401", "0376", "0426"]
filenames = glob("wrfout_d01_2018*", "/home/cfranken/data/WRF-CO2/")

# Dictionary to hold the filenames grouped by timestamp
grouped_files = Dict{String, Vector{String}}()

# Function to parse the filename and extract the timestamp and tile digits
function parse_filename(filename)
    # Match the timestamp and tile part from the filename
    m = match(r"wrfout_d01_(\d{4}-\d{2}-\d{2}_\d{2}_\d{2}_\d{2})_(\d{4})", filename)
    if m !== nothing
        timestamp = m.captures[1]
        tile_digits = m.captures[2]
        return (timestamp, tile_digits)
    else
        return nothing
    end
end

# Dictionary to hold the filenames grouped by timestamp
grouped_files = Dict{String, Vector{String}}()

# Loop through filenames and group them by timestamp if the tile matches
for filename in filenames
    parsed = parse_filename(filename)
    if parsed !== nothing
        timestamp, tile_digits = parsed
        # If the tile digits are in the provided list, add to the group
        if tile_digits in tile_digits_to_include
            if haskey(grouped_files, timestamp)
                push!(grouped_files[timestamp], filename)
            else
                grouped_files[timestamp] = [filename]
            end
        end
    end
end



# Filter out timestamps where not all tiles are present
complete_timesteps = Dict{String, Vector{String}}()
for (timestamp, files) in grouped_files
    # Check if all specified tile digits are present for the timestamp
    tiles_present = [parse_filename(f)[2] for f in files]
    @show tiles_present
    if all(tile -> tile in tiles_present, tile_digits_to_include)
        complete_timesteps[timestamp] = files
    end
end

# Convert the timestamp strings into DateTime objects and sort
sorted_timesteps = sort(collect(keys(complete_timesteps)), by = t -> DateTime(t, "yyyy-mm-dd_HH_MM_SS"))

# Display the result sorted by timestamp
for timestamp in sorted_timesteps
    println("Timestamp: $timestamp")
    println("Files: ", complete_timesteps[timestamp])
    outname = replace(complete_timesteps[timestamp][1], "WRF-CO2/" => "WRF-CO2/plots/")
    outname = outname * ".png"
    @show outname
    extractXCO2(complete_timesteps[timestamp], outname, timestamp)

end




fig = Figure(title="Test")

ga = GeoAxis(fig[1, 1];  limits=((-118.8, -117.5), (33.3, 34.8))) #dest = "+proj=merc",
#lines!(ga, GeoMakie.coastlines(); strokewidth = 0.7, color=:gold, rasterize = 5)
xco2 = GeoMakie.heatmap!(ga,  lons,lats, xco2_interp')
lines!(ga, coasts; color= :grey, xautolimits=false, yautolimits=false, rasterize=5)
# Add a star marker for Pasadena
scatter!(ga,[lon_pasadena], [lat_pasadena], markershape = :star5, markersize = 10, color = :red)
Colorbar(fig[1, 2], xco2, label="WRF-Chem XCO₂ (ppm)", height = Relative(0.85))
#poly!(ga, counties; color= :grey, xautolimits=false, yautolimits=false, rasterize=5)
#poly!(ga, GeoMakie.to_multipoly(coasts.geometry); strokewidth = 0.7, color=:gold, rasterize = 5)
#poly!(ga, GeoMakie.to_multipoly(lakes.geometry); strokewidth = 0.7, color=:blue, rasterize = 5,  xautolimits=false, yautolimits=false)
fig