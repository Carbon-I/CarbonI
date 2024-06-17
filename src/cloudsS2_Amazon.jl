using CSV, DataFrames, Glob, Statistics,  StatsBase
using Makie, CairoMakie

function extract_month(date_string::String)
    # Extracting the month part from the string
    month_str = date_string[5:6]
    # Converting the month string to integer
    parse(Int, month_str)
end

# Function to check if the month is in the range July to September
is_in_dry(row) = row.Month in 7:9
is_in_wet(row) = row.Month in [1,2,3,4,5,12]

x = [30,90,  150, 200, 300, 400, 600, 800, 1000, 1250, 1500, 1750, 2000,3000, 5000, 7500]#, 3000, 5000, 7500, 10000];
datasets = []
for res in (30,90, 150, 200, 300, 400, 600, 800, 1000, 1250, 1500, 1750, 2000,3000, 5000, 7500)#, 3000, 5000, 7500, 10000)
    files = Glob.glob("AmazonBox2018-01-01_2022-01-01_*_"* string(res)*".0.csv","/home/cfranken/shared/GE_BOX2/")
    df = CSV.read(files, DataFrame; source = "source")
    df[!,:Month] = [extract_month(date_string) for date_string in df.var"system:index"]

    push!(datasets,df)
end

med_dry = [quantile(filter(is_in_dry,d).cloud_free_fraction1,[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]) for d in datasets]
med_wet = [quantile(filter(is_in_wet,d).cloud_free_fraction1,[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]) for d in datasets]

cf_mean_dry = [mean(filter(is_in_dry,d).cloud_free_fraction1) for d in datasets]
cf_mean_wet = [mean(filter(is_in_wet,d).cloud_free_fraction1) for d in datasets]

percentiles_wet = mapreduce(permutedims, vcat, med_wet)
percentiles_dry = mapreduce(permutedims, vcat, med_dry)


bins = [0,0.00001,  0.001,0.002,0.005, 0.01,  0.05, 0.1, 0.2, 0.2,0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
bins = [0, 0.001,0.002,0.005, 0.01, 0.03, 0.05, 0.1,  0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
bins = [0; 10.0.^(-3.0:0.1:0)]
hist1_dry = [fit(Histogram, filter(is_in_dry,d).cloud_free_fraction1, bins).weights for d in datasets];
hist5_dry = [fit(Histogram, filter(is_in_dry,d).cloud_free_fraction5, bins).weights for d in datasets];

hist1_wet = [fit(Histogram, filter(is_in_wet,d).cloud_free_fraction1, bins).weights for d in datasets];
hist5_wet = [fit(Histogram, filter(is_in_wet,d).cloud_free_fraction5, bins).weights for d in datasets];

#bla = [fit(Histogram, log10.(d.cloud_free_fraction1), -15:0.1:0).weights for d in datasets]
binsDry1 = mapreduce(permutedims, vcat, hist1_dry)
binsDry5 = mapreduce(permutedims, vcat, hist5_dry)

binsWet1 = mapreduce(permutedims, vcat, hist1_wet)
binsWet5 = mapreduce(permutedims, vcat, hist5_wet)






