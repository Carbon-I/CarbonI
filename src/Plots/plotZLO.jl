using NCDatasets, Statistics, Plots, Glob
sif = Dataset("data/OCO2_LtSIF_180825_B11100Ar_230602034736s.nc4");
off = sif.group["Offset"];
off_757 = off["SIF_Relative_Median_757nm"][:];
sig = off["signal_histogram_bins"][:];


# Define the directory containing the .nc4 files
dir_path = "/home/cfranken/data/SIF_lite/oco2.gesdisc.eosdis.nasa.gov/data/OCO2_DATA/OCO2_L2_Lite_SIF.11r/2016/"

# Get a list of all .nc4 files in the folder, sorted by name
nc4_files = sort(Glob.glob("*.nc4", dir_path))

# Initialize an empty matrix to store the data
combined_matrix_757 = []
combined_matrix_771 = []

# Loop over each file and extract the vector from the specified group and variable
for file in nc4_files
    # Open the .nc4 file
    dataset = Dataset(joinpath(dir_path, file))
    
    # Access the variable in the Offset group
    off_757 = dataset.group["Offset"]["SIF_Relative_Median_757nm"][:]
    off_771 = dataset.group["Offset"]["SIF_Relative_Median_771nm"][:]
    
    # Append the vector as a new row in the growing matrix
    if length(combined_matrix_757) == 0
        combined_matrix_757 = (mean(off_757[2,:,:], dims=1)[1,:])'
        combined_matrix_771 = (mean(off_771[2,:,:], dims=1)[1,:])'
    else
        combined_matrix_757 = vcat(combined_matrix_757, (mean(off_757[2,:,:], dims=1)[1,:])')
        combined_matrix_771 = vcat(combined_matrix_771, (mean(off_771[2,:,:], dims=1)[1,:])')
    end
    
    # Close the dataset after reading
    close(dataset)
end