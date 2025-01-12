using XLSX
using DataFrames
using Plots

include(joinpath(@__DIR__, "Plots", "CI_colors.jl"))

# Load the Excel file
file_path = "data/edgar/EDGAR_CH4_1970_2023.xlsx"
file_path_co2  = "data/edgar/IEA_EDGAR_CO2_1970_2023.xlsx"
raw_data = DataFrame(XLSX.readtable(file_path, "TOTALS BY COUNTRY"))
raw_data_co2 = DataFrame(XLSX.readtable(file_path_co2, "TOTALS BY COUNTRY"))


# Rename columns for clarity
column_names = ["Region", "Country", "Country_Code", "Name", "Substance"]
years = [string("Y_", year) for year in 1970:2023]
rename!(raw_data, vcat(column_names, years))
rename!(raw_data_co2, vcat(column_names, years))

# Filter for specific regions
regions_of_interest = r"(India|China|United States|European)"
filtered_data = filter(row -> occursin(regions_of_interest, string(row.Name)), raw_data)
filtered_data_co2 = filter(row -> occursin(regions_of_interest, string(row.Name)), raw_data_co2)

# Calculate Europe totals
europe_keywords = r"(Europe|EU|European)"
europe_data = filter(row -> occursin(europe_keywords, string(row.Country)), raw_data)
europe_data_co2 = filter(row -> occursin(europe_keywords, string(row.Country)), raw_data_co2)
europe_emissions = convert.(Float64,Matrix(select(europe_data, names(europe_data)[6:end])))
europe_emissions_co2 = convert.(Float64,Matrix(select(europe_data_co2, names(europe_data_co2)[6:end])))
europe_emissions_total = sum(europe_emissions, dims=1)[1,:]
europe_emissions_total_co2 = sum(europe_emissions_co2, dims=1)[1,:]

# Compute global totals and Rest of the World
global_emissions = convert.(Float64,Matrix(select(raw_data, names(raw_data)[6:end]) ))
a = Matrix(select(raw_data_co2, names(raw_data_co2)[6:end]))
a[ismissing.(a)] .= 0
global_emissions_co2 = convert.(Float64, a)
global_emissions_total = sum(global_emissions, dims=1)[1,:]
global_emissions_total_co2 = sum(global_emissions_co2, dims=1)[1,:]


# Combine emissions for plotting
china_emissions = convert.(Float64,Matrix(filtered_data[filtered_data.Name .== "China", 6:end]))[1,:]
china_emissions_co2 = convert.(Float64,Matrix(filtered_data_co2[filtered_data_co2.Name .== "China", 6:end]))[1,:]
india_emissions = convert.(Float64,Matrix(filtered_data[filtered_data.Name .== "India", 6:end]))[1,:]
india_emissions_co2 = convert.(Float64,Matrix(filtered_data_co2[filtered_data_co2.Name .== "India", 6:end]))[1,:]   
us_emissions = convert.(Float64,Matrix(filtered_data[filtered_data.Name .== "United States", 6:end]))[1,:]
us_emissions_co2 = convert.(Float64,Matrix(filtered_data_co2[filtered_data_co2.Name .== "United States", 6:end]))[1,:]
rest_of_world_emissions = global_emissions_total .- (china_emissions .+ us_emissions .+ europe_emissions_total)
rest_of_world_emissions_co2 = global_emissions_total_co2 .- (china_emissions_co2 .+ us_emissions_co2 .+ europe_emissions_total_co2)

# Example data
years = 1970:2023
china =china_emissions/1e3

india = india_emissions/1e3
us = us_emissions/1e3
europe = europe_emissions_total/1e3
rest_of_world = rest_of_world_emissions/1e3

# Compute cumulative emissions for stacking
r = 44/12*1e6
cumulative_rest_of_world = rest_of_world
cumulative_rest_of_world_co2 = rest_of_world_emissions_co2/r
cumulative_us = cumulative_rest_of_world + us
cumulative_us_co2 = cumulative_rest_of_world_co2 + us_emissions_co2/r
cumulative_europe = cumulative_us + europe
cumulative_europe_co2 = cumulative_us_co2 + europe_emissions_total_co2/r
cumulative_china = china + cumulative_europe
cumulative_china_co2 = china_emissions_co2/r + cumulative_europe_co2
cumulative_india = cumulative_china + india
cumulative_india_co2 = cumulative_china_co2 + india_emissions_co2/r



years = 1970:2023
r = 44/12
# Create a figure
fig = Figure(resolution = (800, 700))
ax = Axis(fig[1, 1], title = "Cumulative Emissions Over Time", xlabel = "Year", ylabel = "CH₄ Emissions (Tg CH₄/yr)")
lw = 1
# Plot each layer using `fill_between!`
fill_between!(ax, years, zeros(length(years)), cumulative_rest_of_world, color = (CarbonI_colors[1],0.8), label = "Rest of World", alpha = 0.7)
lines!(ax, years, cumulative_rest_of_world, color = :black, linewidth = lw)
fill_between!(ax, years, cumulative_rest_of_world, cumulative_us, color = (CarbonI_colors[10],0.8), label = "United States", alpha = 0.7)
lines!(ax, years, cumulative_us, color = :black, linewidth = lw)
fill_between!(ax, years, cumulative_us, cumulative_europe, color = (CarbonI_colors[13],0.8), label = "Europe", alpha = 0.7)
lines!(ax, years, cumulative_europe, color = :black, linewidth = lw)
fill_between!(ax, years, cumulative_europe, cumulative_china, color = (CarbonI_colors[4],0.8), label = "China", alpha = 0.7)
lines!(ax, years, cumulative_china, color = :black, linewidth = lw)
fill_between!(ax, years, cumulative_china, cumulative_india, color = (CarbonI_colors[7],0.8), label = "India", alpha = 0.7)
lines!(ax, years, cumulative_india, color = :black, linewidth = lw)


CairoMakie.xlims!(ax, 1980,2023)
CairoMakie.ylims!(ax, 0,400)
ax2 = Axis(fig[2, 1],  xlabel = "Year", ylabel = "CO₂ Emissions (GtC/yr)")
# Plot each layer using `fill_between!`
fill_between!(ax2, years, zeros(length(years)), cumulative_rest_of_world_co2, color = (CarbonI_colors[1],0.8), label = "Rest of World", alpha = 0.7)
lines!(ax2, years, cumulative_rest_of_world_co2, color = :black, linewidth = lw)
fill_between!(ax2, years, cumulative_rest_of_world_co2, cumulative_us_co2, color = (CarbonI_colors[10],0.8), label = "United States", alpha = 0.7)
lines!(ax2, years, cumulative_us_co2, color = :black, linewidth = lw)
fill_between!(ax2, years, cumulative_us_co2, cumulative_europe_co2, color = (CarbonI_colors[13],0.8), label = "Europe", alpha = 0.7)
lines!(ax2, years, cumulative_europe_co2, color = :black, linewidth = lw)
fill_between!(ax2, years, cumulative_europe_co2, cumulative_china_co2, color = (CarbonI_colors[4],0.8), label = "China", alpha = 0.7)
lines!(ax2, years, cumulative_china_co2, color = :black, linewidth = lw)
fill_between!(ax2, years, cumulative_china_co2, cumulative_india_co2, color = (CarbonI_colors[7],0.8), label = "India", alpha = 0.7)
lines!(ax2, years, cumulative_india_co2, color = :black, linewidth = lw)
axislegend(ax2, position=:lt, labelsize =12)
#axislegend(ax2, position=:lt, labelsize =12)
CairoMakie.xlims!(ax2, 1980,2023)
CairoMakie.ylims!(ax2, 0,12)
hidexdecorations!(ax, grid=false)

ax3 = Axis(fig[1, 2]; autolimitaspect=1,title = "Fractional Contributions (and Trend)",)
CairoMakie.pie!(ax3,[rest_of_world[end], us[end], europe[end], china[end], india[end]], colors=[(CarbonI_colors[1],0.8),(CarbonI_colors[10],0.8),(CarbonI_colors[13],0.8),(CarbonI_colors[4],0.8),(CarbonI_colors[7],0.8)])
# Display the plot
display(fig)