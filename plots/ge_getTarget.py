import ee
import time
import xarray as xr
import rioxarray

# Authenticate and initialize Google Earth Engine
ee.Authenticate()  # Authenticate your account
ee.Initialize()    # Initialize the API

# Step 1: Define the Area of Interest (AOI)
# Define the center of the AOI (example: lat=40, lon=-120)
center = ee.Geometry.Point([-120, 40])

# Create a rectangle 100x100 km (50km buffer from the center in each direction)
aoi = center.buffer(50000).bounds()  # 50km buffer in meters from the center

# Step 2: Filter Sentinel-2 Image Collection and Select Band 12
# Sentinel-2 surface reflectance image collection
collection = ee.ImageCollection('COPERNICUS/S2_SR')

# Filter by date range (adjust the dates as needed)
filtered = collection.filterDate('2023-01-01', '2023-01-31')\
                     .filterBounds(aoi)\
                     .select('B12')  # Selecting Band 12 (SWIR)

# Get the median image over the time range and clip to AOI
image = filtered.median().clip(aoi)

# Step 3: Export Image as GeoTIFF to Google Drive
# Export the image to GeoTIFF
export_task = ee.batch.Export.image.toDrive(
    image=image,
    description='Sentinel2_Band12_AOI',
    scale=10,  # Sentinel-2 resolution is 10m for Band 12
    region=aoi.getInfo()['coordinates'],  # Export the AOI
    fileFormat='GeoTIFF',  # Export format
    crs='EPSG:4326',  # Set the coordinate reference system (WGS84)
    maxPixels=1e9  # Adjust depending on the size of your image
)

# Start the export task
export_task.start()

# Step 4: Monitor the Export Task Status
print("Export started, monitoring task status...")

while export_task.status()['state'] in ['RUNNING', 'READY']:
    print(f"Task status: {export_task.status()['state']}")
    time.sleep(30)  # Wait 30 seconds between checks

print(f"Task finished with state: {export_task.status()['state']}")

# Step 5: Convert GeoTIFF to NetCDF (after export completion)
# Load the exported GeoTIFF (replace 'Sentinel2_Band12_AOI.tif' with the path to the downloaded file)
geo_tiff_path = 'Sentinel2_Band12_AOI.tif'  # Adjust this to your actual file path

# Open the GeoTIFF with rioxarray
raster = rioxarray.open_rasterio(geo_tiff_path)

# Convert to NetCDF and save
netcdf_path = 'Sentinel2_Band12_AOI.nc'
raster.to_netcdf(netcdf_path)

print(f"GeoTIFF converted to NetCDF: {netcdf_path}")

