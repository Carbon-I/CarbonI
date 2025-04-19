using Statistics, vSmartMOM, InstrumentOperator, CarbonI, ImageFiltering, Interpolations, Statistics, CUDA, CairoMakie
using DelimitedFiles


Δwl = 0.005
wl = 2035:Δwl:2385

# Define an instrument:
cM, wl_ci = CarbonI.create_carbonI_conv_matrix_cbe(wl)
include("src/Plots/CI_colorsNew.jl")

# Load tropical forest albedo
scenario = CarbonI.stressing_scenario()
# Load the soil albedo data
soil = readdlm("data/soil.aridisol.camborthid.none.all.89p1772.jhu.becknic.spectrum.txt", skipstart=1751)
wo = findall(2.45 .> soil[:,1].> 1.9)
soil = soil[wo,:]

grid = 2000:10:2400
grid_tr = 1980:30:2400
spl = CubicSplineInterpolation(grid_tr, scenario.surface_albedo(grid_tr))
soil_spline = LinearInterpolation(1000*reverse(soil[:,1]), reverse(soil[:,2]))
soil_smooth = CubicSplineInterpolation(grid, soil_spline(grid))

function plotAlbedos()
    f = Figure(resolution=(500,400), title="Spectrally Rsolved Albedos", fontsize=12)
    ax1 = Axis(f[1,1], yminorgridvisible = true,  xlabel="Wavelength (μm)",ylabel="Reflectance (%)", backgroundcolor=:transparent)
    lines!(ax1, wl_ci, scenario.surface_albedo(wl_ci)*100 , label="Aridisol", color=CarbonI_colors[3], linewidth=3)
    lines!(ax1, 1000*soil[:,1], soil[:,2] , label="Aridisol", color=CarbonI_colors[5], linewidth=3)
    f
end

function plotAlbedos()
    f = Figure(resolution=(500,400), title="Spectrally Resolved Albedos", fontsize=12)

    # Primary axis (left Y)
    ax1 = Axis(f[1, 1];
        xlabel="Wavelength (μm)",
        ylabel="Tropical Forest Reflectance (%)",
        yminorgridvisible=true,
        #backgroundcolor=:transparent
    )

    # Twin axis (shares X, right Y)
    ax2 = Axis(f[1, 1];
        #sharex = ax1,                # link X to ax1
        yaxisposition = :right,      # put its Y on the right
        ylabel = "Aridisol Soil Reflectance (%)",
        yminorgridvisible=true,
        #xticksvisible = false,       # hide its X‐axis ticks
        #xticklabelvisible = false,
        #backgroundcolor = :transparent
    )

    # plot on ax1
    lines!(ax1,
        wl_ci,
        spl(wl_ci) .* 100,
        label = "Forest Albedo",
        color = CarbonI_colors[3],
        linewidth = 3,
    )
    lines!(ax2,
    wl_ci,
        soil_smooth(wl_ci),
        label = "Aridisol Soil Albedo",
        color = CarbonI_colors[1],
        linewidth = 3,
    )

    # separate legends so they don’t overlap
    axislegend(ax1, position = :lb)
    axislegend(ax2, position = :rt)

    linkxaxes!(ax1, ax2)
    CairoMakie.xlims!(ax1, 2030, 2380)
    CairoMakie.xlims!(ax2, 2030, 2380)
    CairoMakie.ylims!(ax1, 2.5, 6)
    CairoMakie.ylims!(ax2, 25, 60)
   # xlims!(ax2, 2000, 2400)
    return f
end

f = plotAlbedos()
save("plots/ProxyPaper/AlbedoSpectra.pdf", f)