using NCDatasets, GeoMakie, CairoMakie, Statistics, MAT, Kronecker, Interpolations
# Read bottom-up emissions (from Yi):
bu_ch4 = Dataset("data/GCP_wetland_methane.nc")
bu_mat = matread("data/bottomup_ch4.mat")
lat_mat = bu_mat["lat"];
lon_mat = bu_mat["lon"];
ch4_wetchimp_unc = bu_mat["ch4_wetchimp_unc"]'[:,end:-1:1];
ch4_wetchimp = bu_mat["ch4_wetchimp_mean"]'[:,end:-1:1];
ch4_gcp_unc = bu_mat["ch4_gcp_unc"]'[:,end:-1:1];
ch4_gcp = bu_mat["ch4_gcp_mean"]'[:,end:-1:1];
ch4_wetcharts_unc = bu_mat["ch4_wetcharts_unc"]'[:,end:-1:1];
ch4_wetcharts = bu_mat["ch4_wetcharts_mean"]'[:,end:-1:1];

# Read Data yield (from Yinon):
data_yield_1     = Dataset("/home/ymbaron/data/cloud_effect/results/01_calc_OCO_yield/OCO2_L2_Lite_FP_11r_2020_binned.nc")
data_yield_raw_1 = Dataset("/home/ymbaron/data/cloud_effect/results/01_calc_OCO_yield/OCO2_L1B_Science_11r_2020_binned.nc")
data_yield_2     = Dataset("/home/ymbaron/data/cloud_effect/results/01_calc_OCO_yield/OCO2_L2_Lite_FP_11r_2019_binned.nc")
data_yield_raw_2 = Dataset("/home/ymbaron/data/cloud_effect/results/01_calc_OCO_yield/OCO2_L1B_Science_11r_2019_binned.nc")
data_yield_3     = Dataset("/home/ymbaron/data/cloud_effect/results/01_calc_OCO_yield/OCO2_L2_Lite_FP_11r_2021_binned.nc")
data_yield_raw_3 = Dataset("/home/ymbaron/data/cloud_effect/results/01_calc_OCO_yield/OCO2_L1B_Science_11r_2021_binned.nc")

data_net = sum(data_yield_1["__xarray_dataarray_variable__"][:], dims=3)[:,:,1] + sum(data_yield_2["__xarray_dataarray_variable__"][:], dims=3)[:,:,1] + sum(data_yield_3["__xarray_dataarray_variable__"][:], dims=3)[:,:,1];
data_raw = sum(data_yield_raw_1["__xarray_dataarray_variable__"][:], dims=3)[:,:,1] + sum(data_yield_raw_2["__xarray_dataarray_variable__"][:], dims=3)[:,:,1] + sum(data_yield_raw_3["__xarray_dataarray_variable__"][:], dims=3)[:,:,1];



# Compute Data Yield as a percentage
data_yield_all = 100 .* (data_net ./ data_raw);
data_yield_0_5deg = kron(data_yield_all, ones(2,2))
# Scatter plot
fig, ax, sc = scatter(data_yield[:], bu_ch4["STD"][:][:], markersize=3, color=:blue)

# Labels
ax.xlabel = "Data Yield (%)"
ax.ylabel = "CH₄ Emissions (STD)"

# Display the figure
fig

cmap = cgrad(:RdYlBu_11, 10, rev=true, categorical = true)
cmap2 = cgrad(:RdYlBu_11, 10,  categorical = true)
#cmap = :RdYlBu_11
ch4_gcp_unc[ch4_gcp_unc .< 0.01] .=NaN

ch4_unc_mean = (ch4_gcp_unc .+ ch4_wetcharts_unc .+ ch4_wetchimp_unc)./3;
data_yield_0_5deg[isnan.(ch4_unc_mean)] .= NaN

levels = [0.0, 0.5, 1, 2, 5, 10, 15, 20, 30, 50, 70]
cScale = LinearInterpolation(levels, collect(1:length(levels)), extrapolation_bc=NaN)

function plot_data_yield()
    fig = Figure(resolution=(800,400))
    ga = GeoAxis(title="OCO-2 Data Yield (%)",
        fig[1, 1], limits=((-180.0, 180.0), (-60, 80)),
        xlabelvisible=false, ylabelvisible=false,  # Hide lat/lon labels
        xticklabelsvisible=false, yticklabelsvisible=false, 
        xticksvisible=false, yticksvisible=false); 

    #sp = GeoMakie.heatmap!(ga, lon, lat, (sum(N, dims=1)[1,:,:]); colorrange=(10^1,10^5), colormap=cmap)
    sp = GeoMakie.heatmap!(ga, lon_mat, lat_mat, cScale.(data_yield_0_5deg); colormap=cmap2, colorrange=(1,11));
    cb = Colorbar(fig[2, 1], sp;  width = Relative(0.85), vertical=false)
    cb.ticks[] = (collect(1:length(levels)), string.(levels))
    lines!(ga, GeoMakie.coastlines(), color=:gray, linewidth=2) 

    fig
end


function plot_ch4_unc()
    levels2 = [0.0, 3, 6, 10, 20, 30, 50, 75, 100, 150, 200]
    cScale2 = LinearInterpolation(levels2, collect(1:length(levels2)), extrapolation_bc=NaN)
    fig = Figure(resolution=(800,400))
    ga = GeoAxis(title="Mean CH₄ Emission Uncertainties (mg/m²/day)",
        fig[1, 1], limits=((-180.0, 180.0), (-60, 80)),
        xlabelvisible=false, ylabelvisible=false,  # Hide lat/lon labels
        xticklabelsvisible=false, yticklabelsvisible=false, 
        xticksvisible=false, yticksvisible=false); 

    #sp = GeoMakie.heatmap!(ga, lon, lat, (sum(N, dims=1)[1,:,:]); colorrange=(10^1,10^5), colormap=cmap)
    sp = GeoMakie.heatmap!(ga, lon_mat, lat_mat, cScale2.(1000/365*ch4_unc_mean); colormap=cmap, colorrange=(1,11));
    cb = Colorbar(fig[2, 1], sp;  width = Relative(0.85), vertical=false)
    cb.ticks[] = (collect(1:length(levels2)), string.(levels2))
    lines!(ga, GeoMakie.coastlines(), color=:gray, linewidth=2) 

    fig
end
fig = plot_data_yield()
CairoMakie.save("plots/oco2_yield.pdf", fig)
fig = plot_ch4_unc()
CairoMakie.save("plots/ch4_uncertainties.pdf", fig)

fig2 = with_theme(plot_data_yield, theme_black())
CairoMakie.save("plots/oco2_yield_dark.pdf", fig2)

fig2 = with_theme(plot_ch4_unc, theme_black())
CairoMakie.save("plots/ch4_uncertainties_dark.pdf", fig2)