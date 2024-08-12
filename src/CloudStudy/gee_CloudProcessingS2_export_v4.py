#!/usr/bin/python3

import argparse
import geemap
import ee
import datetime
import numpy as np
import h5py
from datetime import timedelta
from dateutil.relativedelta import relativedelta
import pandas as pd
import time
import progressbar
import sys
import os

# Initialize the Earth Engine
ee.Initialize()

# Define cloud mask used:
def cloud_mask(image):
        cloud = image.select(['cs_cdf']).lt(0.65)
        return cloud

# Get Image ID
def get_image_id(image):
    return ee.String(image.id())

# Count number of pixels in basin (as prefilter)
def count_collection_pixels(images, basin, name):
    def count_pixels(image):
        # Mask the image with the domain and count the number of pixels
        count = image.clip(basin).reduceRegion(
            reducer=ee.Reducer.count(),
            geometry=basin,
            maxPixels=1e10
        ).get(name)  # Replace 'B1' with a relevant band name from your dataset
        return image.set('domain_pixel_count', count)
    return images.map(count_pixels)

# Reduce S2 data to Landsat resolution as starting point
def reduceToLandsat(images, bands):
    #resampleScale = fp;
    resampleScale = 30;
    def reduceFunction(image):
        return image.reduceResolution(reducer=ee.Reducer.mean(),maxPixels=64000).reproject(crs=image.projection().crs(),scale=resampleScale)

    return images.map(reduceFunction)

# Convolve with square kernel:
def reduceAll(images, bands, fp, basin):
    kernel=ee.Kernel.square(fp/2, units='meters')
    def reduceFunction(image):
        return image.clip(basin).reduceNeighborhood(reducer=ee.Reducer.sum(),optimization='boxcar',kernel=kernel).rename(bands)#.reproject(crs=image.projection().crs(),scale=fp).rename(bands)#.reproject(image.projection().atScale(resampleScale))#.sample(scale=fp)
        #return image.clip(basin).reduceResolution(reducer=ee.Reducer.mean(),maxPixels=64000).reproject(crs=image.projection().crs(),scale=fp)
        #retucloud_mask_imagern image.clip(basin).convolve(kernel=kernel).rename(bands)#.reproject(image.projection().atScale(resampleScale))#.sample(scale=fp)
    return images.map(reduceFunction)

def fp_loops(cloud_mask_images, footprint_sizes, max_cloud_percentages, basin):
    def calc_cloud_free_frac_loop(cloud_mask_image):
        result = ee.Feature(None,{})#ee.List([])
        for footprint in footprint_sizes:
            #print(footprint)
            #try:
            # reduce the resolution of the cloud mask to the footprint scale and calculate the fraction of cloud free pixels in the footprint
            cloud_fraction_map = cloud_mask_image.reduceResolution(reducer=ee.Reducer.mean(),
                                    maxPixels=65535).reproject(
                                    crs=cloud_mask_image.projection().crs(),
                                    scale=footprint)
            cloud_mask_image = cloud_fraction_map
            # calculate the total number of pixels in the footprint
            total_pixels = cloud_fraction_map.reduceRegion(
                        reducer=ee.Reducer.count(),
                        geometry=basin,
                        maxPixels=1e13
            )
            
            # Loop through fractions
            for cloud_threshold in max_cloud_percentages:
                # calculate the total number of pixels with a lower fraction of clouds than the threshold
                cloud_free_pixels = cloud_fraction_map.lt(cloud_threshold).reduceRegion(
                    reducer=ee.Reducer.sum(),
                    geometry=basin,
                    maxPixels=1e10
                )
                # calculate the fraction of cloud free pixels in the footprint
                cloud_free_fraction = ee.Number(cloud_free_pixels.get('cs_cdf')).divide(ee.Number(total_pixels.get('cs_cdf')))
                result = result.set(f'{footprint}_{int(cloud_threshold*100)}',cloud_free_fraction)
            #except:
            #    print("caught error")
                
        return result#image.set(f'frac',result)
    return ee.FeatureCollection(cloud_mask_images.map(calc_cloud_free_frac_loop))

# Process cloud fraction info and store in Image itself:s
def processCloudImages(cloud_images, basin, name):
    def calc_cloud_free_frac(image):
        cloud_free_pixels1 = image.lt(0.01).reduceRegion(
                                    reducer=ee.Reducer.sum(),
                                    geometry=basin,
                                    maxPixels=1e13
                                ).get(name)
        #cloud_free_pixels10 = image.lt(0.1).reduceRegion(
        #                            reducer=ee.Reducer.sum(),
        #                            geometry=basin,
        #                            maxPixels=1e10
        #                        ).get(name)
        cloud_free_pixels5 = image.lt(0.002).reduceRegion(
                                    reducer=ee.Reducer.sum(),
                                    geometry=basin,
                                    maxPixels=1e13
                                ).get(name)
        pixels = image.reduceRegion(
                                    reducer=ee.Reducer.count(),
                                    geometry=basin,
                                    maxPixels=1e13
                                ).get(name)
        
        #cloud_free_fraction05 = ee.Number(cloud_free_pixels05).divide(ee.Number(pixels))
        cloud_free_fraction1 = ee.Number(cloud_free_pixels1).divide(ee.Number(pixels))
        #cloud_free_fraction10 = ee.Number(cloud_free_pixels10).divide(ee.Number(pixels))
        cloud_free_fraction5 = ee.Number(cloud_free_pixels5).divide(ee.Number(pixels))
        #image = image.set(f'cloud_free_fraction05',cloud_free_fraction05)
        image = image.set(f'cloud_free_fraction1',cloud_free_fraction1)
        #image = image.set(f'cloud_free_fraction10',cloud_free_fraction10)
        image = image.set(f'cloud_free_fraction02',cloud_free_fraction5)
        image = image.set(f'domain_pixel_count',pixels)
        #image = image.set('domain_pixel_count', total_pixels)
        return image
    return cloud_images.map(calc_cloud_free_frac)

def count_ready_tasks():
    """Count the number of tasks in the 'READY' state."""
    tasks = ee.batch.Task.list()
    ready_tasks = [task for task in tasks if task.state == 'READY']
    return len(ready_tasks)

def calculate_number_of_steps(start_date, end_date, step_size):
    """Calculate the number of steps between start and end dates with given step size."""
    total_months = (end_date.year - start_date.year) * 12 + end_date.month - start_date.month
    return (total_months + step_size - 1) // step_size  # Round up to include end date

def processAll(args):

    # Basic stuff
    current_date           = args.start_date
    footprint_sizes        = args.footprints
    max_cloud_percentages  = args.cloud_fractions
    number_of_steps = calculate_number_of_steps(args.start_date, args.end_date, args.step)
    start_lat = np.arange( args.latMin, args.latMax, args.dLat)
    start_lon = np.arange( args.lonMin, args.lonMax, args.dLon)
    
    nAll = number_of_steps * len(start_lat) * len(start_lon)
    iterI = 0
    # Print an initial message before the loop
    mes = "Starting the process... Number of time steps: " + str(number_of_steps) +  " number of latitudes: " + str(len(start_lat)) +  " number of longitudes: " + str(len(start_lon)) +  "All:" + str(nAll)
    print(mes)
    
    bands = ['cs_cdf'];
    num_ready_tasks = 0
    for i in range(number_of_steps):
        current_date      = args.start_date + relativedelta(months=i * args.step)
        current_end_date  = args.start_date + relativedelta(months=(i+1) * args.step)
        # loop through regions
        for iJ in range(0,len(start_lat)):
            for iI in range(0,len(start_lon)):
                iterI +=1
                
                # Start the high-resolution timer
                start_time = time.perf_counter()
                print(start_lat[iJ],start_lon[iI], current_date, iterI/nAll*100, "%")
                filename_check = "/home/cfranken/shared/GE_output_test/" + args.outFile + current_date.strftime("%Y-%m-%d") + "_"  + current_end_date.strftime("%Y-%m-%d") + "_" + str(start_lon[iI]) + "_" + str(start_lat[iJ]) +"_1000.0.csv"
                print(filename_check)
                if os.path.exists(filename_check):
                    print("Skipping file ", filename_check)
                    continue

                # Set region (make sure overlap is in center of the basin, i.e. basin2)
                diffi = min(args.dLat/3.5,0.6)
                basin = ee.Geometry.Rectangle([start_lon[iI], start_lat[iJ], start_lon[iI]+args.dLon, start_lat[iJ]+args.dLat])
                basin2 = ee.Geometry.Rectangle([start_lon[iI]+diffi, start_lat[iJ]+diffi, start_lon[iI]+args.dLon-diffi, start_lat[iJ]+args.dLat-diffi])
                combined_data = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED')
                
                # Filter area, time, and land only:
                filtered_data = combined_data.filterBounds(basin2).select(bands).filterDate(current_date.strftime("%Y-%m-%d"), current_end_date.strftime("%Y-%m-%d"))
                pixel_counted_collection = count_collection_pixels(filtered_data, basin, 'cs_cdf')

                # Careful, this depends a bit on the size of the region. this value works well for 1 degree!
                min_pixels = 1.0e7 * np.cos(start_lat[iJ]/180*np.pi)

                # Make sure enough data is in each pixel and limit to a certain total amount:
                filtered_collection = pixel_counted_collection.filterMetadata('domain_pixel_count', 'greater_than', min_pixels)
                
                #print("Number of potential images: ", filtered_collection.size().getInfo())
                #filtered_collection = filtered_collection.limit(12)
                try:
                    
                    # Start Cloud analysis:
                    cloud_mask_images = filtered_collection.map(cloud_mask)

                    # Reduce to 30m as starting point:
                    cloud_mask_images_30 = reduceToLandsat(cloud_mask_images, 'cs_cdf')
                        
                    # Loop until the number of 'READY' tasks is less than 2500
                    while True:
                        #if num_ready_tasks < 100:
                        #    num_ready_tasks += 2
                        #    break

                        # Only check again if previous one was 
                        num_ready_tasks = count_ready_tasks()
                        print(f"Number of 'READY' tasks: {num_ready_tasks}")

                        if num_ready_tasks < 400:
                            break

                        print("Waiting for tasks to start...")
                        time.sleep(10)  # Wait for 20 seconds before checking again
                        
                    # Loop through footprints 
                    for idF, footprint in enumerate(args.footprints):
                        # Don't downsample the 30m run!
                        if footprint > 40:
                            downSample1 = reduceAll(cloud_mask_images_30, bands,footprint, basin);
                        else:
                            downSample1 = cloud_mask_images_30

                            
                        # Loop over cloud thresholds
                        # for idCF, cf in enumerate(args.cloud_fractions):
                        cloud_frac_res = processCloudImages(downSample1,basin,'cs_cdf')
                        filename = args.outFile + current_date.strftime("%Y-%m-%d") + "_"  + current_end_date.strftime("%Y-%m-%d") + "_" + str(start_lon[iI]) + "_" + str(start_lat[iJ]) +"_" +  str(footprint)
                        filename_check = "/home/cfranken/code/gitHub/GE_Tropics/" + current_date.strftime("%Y-%m-%d") + "_"  + current_end_date.strftime("%Y-%m-%d") + "_" + str(start_lon[iI]) + "_" + str(start_lat[iJ]) +"_" +  str(footprint) + ".csv"

                            
                        # Export as batch process
                        task = ee.batch.Export.table.toDrive(
                                    collection=cloud_frac_res,
                                    description=filename,
                                    folder = args.outFolder,
                                    fileFormat='CSV',
                        )
                        ee.Algorithms.If(ee.Number(filtered_collection.size()).gt(0),
                            task.start()
                        )
                    end_time = time.perf_counter()
                    elapsed_time = end_time - start_time
                    print(elapsed_time, "s")
                    sys.stdout.flush()
                    ##message = " Lat:" +  str(start_lat[iJ]) + ", Lon:" + str(start_lon[iI]) +  ", Time:" +  str(current_date) +  ", Image size:" + str(result_df.size) +  ", Time elapsed (s): " + str(elapsed_time)
                    #print(message)
                except Exception as e:
                    print("Error caught")
                    print(f"Error: {str(e)}")
                    sys.stdout.flush()
    
def standalone_main():
    # Use new argument parser:
    parser = argparse.ArgumentParser(description='Run google Earth Engine Landsat cloud statistics. Save information in output HDF5 file')
    parser.add_argument('--footprints', dest='footprints', type=float, nargs='+', default=[30, 60, 120, 200, 400, 800, 1500, 2000, 3000],
                   help='Footprint sizes for the analysis')
    parser.add_argument('--cloud_fractions', dest='cloud_fractions', type=float, nargs='+', default=[0.002,0.01],
                   help='Cloud fraction thresholds for the analysis')
    parser.add_argument('--dLat', dest='dLat', type=float, default=1.0,
                   help='Latitude bin size')
    parser.add_argument('--dLon', dest='dLon', type=float, default=1.0,
                   help='Longitude bin size')
    parser.add_argument('--latMin', dest='latMin', type=float,default=1.0, 
                   help='Minimum latitude')
    parser.add_argument('--latMax', dest='latMax', type=float, default=2.0, 
                   help='Maximum latitude')
    parser.add_argument('--lonMin', dest='lonMin', type=float, default=-55.0, 
                   help='minimum longitude')
    parser.add_argument('--lonMax', dest='lonMax', type=float, default=-54.0, 
                   help='maximum longitude')
    parser.add_argument('--start_date', dest='start_date', type=datetime.date.fromisoformat)
    parser.add_argument('--end_date', dest='end_date', type=datetime.date.fromisoformat)
    
    parser.add_argument('--outFile', dest='outFile', default='Amazon_Test_',
                   help='Start of the output files')
    parser.add_argument('--outFolder', dest='outFolder', default='GE_output_test',
                   help='Folder for storage')
    parser.add_argument('--credentials', dest='credentials', default='~/.config/earthengine/credentials',
                   help='Earth Engine credentials')
    parser.add_argument("--step", dest='step', type=int, help="Time step size in months", default=1)

    args = parser.parse_args()
    print(args.footprints)
    processAll(args)
if __name__ == "__main__":
    standalone_main()
