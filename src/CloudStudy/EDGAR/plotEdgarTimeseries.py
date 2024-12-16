import pandas as pd
import matplotlib.pyplot as plt

# Load the Excel file
file_path = 'data/edgar/EDGAR_CH4_1970_2023.xlsx'

# Load the relevant sheet and skip metadata rows
data_totals_by_country = pd.read_excel(file_path, sheet_name='TOTALS BY COUNTRY', skiprows=8)

# Rename columns for clarity
data_totals_by_country.columns = ['Region', 'Country', 'Country_Code', 'Name', 'Substance'] + \
                                 [f"Y_{year}" for year in range(1970, 2024)]

# Filter for specific regions
regions_of_interest = ['India', 'China', 'United States', 'European Union']
filtered_data = data_totals_by_country[data_totals_by_country['Name'].str.contains('|'.join(regions_of_interest), na=False)]

# Compute totals for Europe
europe_keywords = ['Europe', 'EU', 'European']
europe_data = data_totals_by_country[data_totals_by_country['Name'].str.contains('|'.join(europe_keywords), na=False)]
europe_emissions = europe_data[[col for col in europe_data.columns if col.startswith('Y_')]].sum()

# Compute global emissions and Rest of World
global_emissions = data_totals_by_country[[col for col in data_totals_by_country.columns if col.startswith('Y_')]].sum()
specified_regions = filtered_data[[col for col in filtered_data.columns if col.startswith('Y_')]].sum()
rest_of_world_emissions = global_emissions - (specified_regions + europe_emissions)

# Combine emissions for plotting
regions_totals = pd.DataFrame({
    'Year': range(1970, 2024),
    'China': filtered_data[filtered_data['Name'] == 'China'].iloc[:, 5:].values.flatten(),
    'India': filtered_data[filtered_data['Name'] == 'India'].iloc[:, 5:].values.flatten(),
    'United States': filtered_data[filtered_data['Name'] == 'United States'].iloc[:, 5:].values.flatten(),
    'Europe': europe_emissions.values,
    'Rest of World': rest_of_world_emissions.values
})

# Ensure all values are numeric
for region in ['China', 'India', 'United States', 'Europe', 'Rest of World']:
    regions_totals[region] = pd.to_numeric(regions_totals[region], errors='coerce').fillna(0)

# Prepare data for plotting
years = regions_totals['Year'].values
emissions = [regions_totals[region].values for region in ['China', 'India', 'United States', 'Europe', 'Rest of World']]

# Plot a stacked area chart
fig, ax = plt.subplots(figsize=(12, 8))
ax.stackplot(years, emissions, labels=['China', 'India', 'United States', 'Europe', 'Rest of World'], alpha=0.7)

# Customize the plot
ax.set_title('Sectoral CH4 Emissions Over Time', fontsize=16)
ax.set_xlabel('Year', fontsize=14)
ax.set_ylabel('Emissions (Gg)', fontsize=14)
ax.legend(title='Region', loc='upper left', fontsize=12)
plt.tight_layout()

plt.show()