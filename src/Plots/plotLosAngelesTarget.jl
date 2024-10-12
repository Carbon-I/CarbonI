using CairoMakie
using ArchGDAL
using Colors
using Interpolations

# Step 1: Load the GeoTIFF
function load_geotiff(filepath)
    ArchGDAL.open(filepath) do dataset
        band = ArchGDAL.read(dataset, 1)  # Read the first band (assuming albedo is in band 1)
        geotransform = ArchGDAL.gettransform(dataset)  # Correct method to get the geotransform
        array = Matrix(band)  # Convert the band to a Julia matrix
        return array, geotransform
    end
end

# Step 2: Plot the GeoTIFF in BW
function plot_geotiff_bw(image_array)
    # Normalize image for BW display (0 to 1 range)
    norm_image = (image_array .- minimum(image_array)) ./ (maximum(image_array) - minimum(image_array))

    # Display the image in black and white using CairoMakie
    fig, ax = Figure(), Axis(fig[1, 1], title = "GeoTIFF Image (Los Angeles)")
    heatmap!(ax, norm_image, colormap = :gray, colorrange = (0, 1))
    return fig
end

# Step 3: Extract Albedo Cross-Section at a Given Latitude
function extract_cross_section(image_array, geotransform, target_latitude)
    x_min, x_res, _, y_max, _, y_res = geotransform
    # Calculate the pixel corresponding to the given latitude
    target_row = Int(floor((y_max - target_latitude) / -y_res))  # Adjust for negative y_res

    # Extract albedo values along the row
    albedo_values = image_array[target_row, :]
    longitudes = x_min .+ (0:size(albedo_values, 1) - 1) * x_res

    return longitudes, albedo_values
end

# Step 4: Plot Albedo Cross-Section
function plot_cross_section(longitudes, albedo_values, target_latitude)
    fig, ax = Figure(), Axis(fig[1, 1], title = "Albedo Cross-Section at Latitude: $target_latitude")
    lines!(ax, longitudes, albedo_values, color = :black, label = "Albedo")
    ax.xlabel = "Longitude"
    ax.ylabel = "Albedo"
    return fig
end



# Example usage
filepath = "data/Sentinel2_Band12_AOI_LosAngeles.tif"
target_latitude = 34.05  # Example latitude for Los Angeles

# Load the GeoTIFF
geoarray = GeoArrays.read(filepath)

# Plot the GeoTIFF in BW
fig1 = plot_geotiff_bw(image_array)
display(fig1)

# Extract the albedo cross-section at the given latitude
longitudes, albedo_values = extract_cross_section(image_array, geotransform, target_latitude)

# Plot the albedo cross-section
fig2 = plot_cross_section(longitudes, albedo_values, target_latitude)
display(fig2)
