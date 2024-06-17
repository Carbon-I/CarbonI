import ee
import numpy as np
import pandas as pd
import seaborn as sns
ee.Initialize()

def analyze_image_collection(image_collection, footprint,cloud_threshold=0.01):
    def cloud_mask(image):
        # create a cloud mask for the image
        qa = image.select('QA_PIXEL')
        cloud = qa.bitwiseAnd(1 << 3).Or(qa.bitwiseAnd(1 << 2)).Or(qa.bitwiseAnd(1 << 4))
        return cloud

    def calc_cloud_free_frac(image):
        # get the cloud mask
        cloud_mask_image = cloud_mask(image)
        # reduce the resolution of the cloud mask to the footprint scale and calculate the fraction of cloud free pixels in the footprint
        cloud_fraction_map = cloud_mask_image.reduceResolution(reducer=ee.Reducer.mean(),
                                maxPixels=1e4).reproject(
                                    crs=cloud_mask_image.projection().crs(), 
                                    scale=footprint)
        
        # calculate the total number of pixels with a lower fraction of clouds than the threshold
        cloud_free_pixels = cloud_fraction_map.lt(cloud_threshold).reduceRegion(
                    reducer=ee.Reducer.sum(),
                    geometry=amazon_basin,
                    maxPixels=1e10
                )
        # calculate the total number of pixels in the footprint
        total_pixels = cloud_fraction_map.reduceRegion(
                    reducer=ee.Reducer.count(),
                    geometry=amazon_basin,
                    maxPixels=1e10
                )
        # calculate the fraction of cloud free pixels in the footprint
        cloud_free_fraction = ee.Number(cloud_free_pixels.get('QA_PIXEL')).divide(ee.Number(total_pixels.get('QA_PIXEL')))
        return image.set('frac',cloud_free_fraction)
    
    # calculate the cloud free fraction for each image in the collection
    cloud_free_map = image_collection.map(calc_cloud_free_frac)
    return cloud_free_map


# Define the area of interest: Amazon Basin
amazon_basin = ee.Geometry.Rectangle([-78.0, -12.0, -48.0, 2.0])

# Load Landsat 8 and 9 image collections
landsat8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(amazon_basin)
landsat9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(amazon_basin)
combined_landsat = landsat8.merge(landsat9)

# Filter by date
# filtered_landsat = combined_landsat.filterDate('2023-01-01', '2023-12-31')
filtered_landsat = combined_landsat.filterDate('2023-01-01', '2023-01-02')

# Compute cloud-free likelihood for different footprint sizes and cloud cover thresholds
footprint_sizes = [60, 90, 120, 150, 300, 600, 900, 1200, 2000, 3000]  # Footprint sizes in meters
max_cloud_percentages = [0.01]  # Maximum cloud cover percentages

results = pd.DataFrame(index=np.arange(filtered_landsat.size().getInfo()),columns=footprint_sizes) #create results df

# run analysis over footprint sizes
for size in footprint_sizes:
    for max_cloud_percentage in max_cloud_percentages:
        likelihood = np.array(analyze_image_collection(filtered_landsat, size,max_cloud_percentage).aggregate_array('frac').getInfo())
        results[size] = likelihood
        print(f'Likelihood of less than {max_cloud_percentage*100}% cloud cover for footprint size {size}m: ', likelihood.mean().round(2),'±',likelihood.std().round(2))

# plot results
results.index.name = 'sample'
results.columns.name = 'footprint size (m)'
data = results.stack().reset_index(name='cloud free pixel fraction')
ax = sns.lineplot(data=data, x='footprint size (m)', y='cloud free pixel fraction')
ax.set(ylim=[0,0.5])