using DelimitedFiles
using Interpolations
using Distributions, CairoMakie, ColorSchemes, LaTeXStrings
#using LaTeXStrings
#using Plots
using vSmartMOM
using vSmartMOM.Scattering
using CarbonI

include("src/Plots/CI_colorsNew.jl")
clima_alb = readdlm(CarbonI.albedo_file,',', skipstart=1)

FT = Float64
# Output parameters
wl = 400.:50.:2500.
vμ = 0.1:0.1:1.0
vσ = 1.5:0.25:2.0

dpath="/home/sanghavi/data/Carbon-I/aerosol/hitran_ri/" #ascii/single_files/"
#fname1="sutherland_khanna_biomass.dat"
fname1 = "shettle_1.dat" 
S1 = readdlm(dpath*fname1)

wl0 = S1[6:end, 1] * 1e3 # μm->nm
wo = findall(x -> 380.0<=x<=2500.0, wl0)
wl1               = wl0[wo]
water_real_ri     = S1[6:end, 2][wo]
water_imag_ri     = S1[6:end, 3][wo]
ice_real_ri       = S1[6:end, 4][wo]
ice_imag_ri       = S1[6:end, 5][wo]
# NaCl_real_ri      = S1[6:end, 6][wo]
# NaCl_imag_ri      = S1[6:end, 7][wo]
seasalt_real_ri   = S1[6:end, 8][wo]
seasalt_imag_ri   = S1[6:end, 9][wo]

fname1 = "shettle_2.dat" 
S1 = readdlm(dpath*fname1)

wl0 = S1[6:end, 1] * 1e3 # μm->nm
wo = findall(x -> 380.0<=x<=2500.0, wl0)
wl2                    = wl0[wo]
watersol_real_ri       = S1[6:end, 2][wo]
watersol_imag_ri       = S1[6:end, 3][wo]
NH2SO4_real_ri         = S1[6:end, 4][wo]
NH2SO4_imag_ri         = S1[6:end, 5][wo]
carbonaceous_real_ri   = S1[6:end, 6][wo]
carbonaceous_imag_ri   = S1[6:end, 7][wo]

fname1 = "shettle_3.dat" 
S1 = readdlm(dpath*fname1)

wl0 = S1[6:end, 1] * 1e3 # μm->nm
wo = findall(x -> 380.0<=x<=2500.0, wl0)
wl3                    = wl0[wo]
volcdust_real_ri       = S1[6:end, 2][wo]
volcdust_imag_ri       = S1[6:end, 3][wo]
H2SO4_real_ri          = S1[6:end, 4][wo]
H2SO4_imag_ri          = S1[6:end, 5][wo]

fname1 = "shettle_5.dat" 
S1 = readdlm(dpath*fname1)

wl0 = S1[6:end, 1] * 1e3 # μm->nm
wo = findall(x -> 380.0<=x<=2500.0, wl0)
wl5                    = wl0[wo]
dustlike_real_ri       = S1[6:end, 6][wo]
dustlike_imag_ri       = S1[6:end, 7][wo]

aer_types = ["water", "ice", "seasalt", 
    "watersol", "NH2SO4", "carbonaceous", 
    "volcdust", "H2SO4", "dustlike"]
Naer = length(aer_types) # Number of aerosol types
real_ri_grid = [water_real_ri', ice_real_ri', seasalt_real_ri', watersol_real_ri', NH2SO4_real_ri', carbonaceous_real_ri', volcdust_real_ri', H2SO4_real_ri', dustlike_real_ri'] 
imag_ri_grid = [water_imag_ri', ice_imag_ri', seasalt_imag_ri', watersol_imag_ri', NH2SO4_imag_ri', carbonaceous_imag_ri', volcdust_imag_ri', H2SO4_imag_ri', dustlike_imag_ri']

re_aer_interp = []
im_aer_interp = []

real_ri_wl = zeros(Naer, length(wl))
imag_ri_wl = zeros(Naer, length(wl))

for i=1:Naer
    re_interp_tmp = LinearInterpolation(wl1, real_ri_grid[i,:][1]')
    push!(re_aer_interp, re_interp_tmp)
    im_interp_tmp = LinearInterpolation(wl1, imag_ri_grid[i,:][1]')
    push!(im_aer_interp, im_interp_tmp)
    for j=1:length(wl)
        real_ri_wl[i,j] = re_aer_interp[i](wl[j])
        imag_ri_wl[i,j] = im_aer_interp[i](wl[j])
    end
end



#======================================================================#
# Simpler size distribution and n_i scheme                             #
#======================================================================#
vμ = [0.08,0.13,0.2] # microns
vσ = [1.3, 1.6]
Niaer=[6,8]
nquad_radius = 2500
r_max = 60.0 #μm
#for ctr=1:9
    #make_mie_model()

extXS   = zeros(length(Niaer), length(vμ), length(vσ), length(wl))
scattXS = zeros(length(Niaer), length(vμ), length(vσ), length(wl))
extC    = zeros(length(Niaer), length(vμ), length(vσ), length(wl))
scattC  = zeros(length(Niaer), length(vμ), length(vσ), length(wl))

for iμ = 1:length(vμ)
    μ = vμ[iμ]
    for iσ = 1:length(vσ)
        σ = vσ[iσ]
        size_distribution = LogNormal(log(μ), log(σ))
        for i=1:length(wl)
            for ictr=1:length(Niaer)#iaer = 1:Naer
                iaer = Niaer[ictr]
                #real_ri = real_ri_wl[iaer, i]
                #imag_ri = imag_ri_wl[iaer, i]
                real_ri = real_ri_wl[iaer, i]
                imag_ri = 0.001 #0.001 #imag_ri_wl[iaer, i]
                aerosol=vSmartMOM.Scattering.Aerosol(size_distribution, real_ri, imag_ri)
                extXS[ictr, iμ, iσ, i], 
                scattXS[ictr, iμ, iσ, i], 
                extC[ictr, iμ, iσ, i], 
                scattC[ictr, iμ, iσ, i] = 
                    vSmartMOM.Scattering.compute_aerosol_XS(aerosol, wl[i]/1e3, r_max, nquad_radius)
            end
        end
    end
end


CS = ColorSchemes.seaborn_colorblind
include(joinpath(@__DIR__, "Plots", "CI_colors.jl"))


aer_size = 0.01:0.001:3
n = 8 #normalization wavelength

function plotAeroSizes()
    f = Figure(resolution=(800,400), title="Aerosol Microphysics", fontsize=18, backgroundcolor = :transparent,fonts = (; regular = "Helvetica Condensed Light", bold="Helvetica Condensed Bold"))
    ax1 = Axis(f[1,2], yminorgridvisible = true,  backgroundcolor = :transparent,yminorticks = IntervalsBetween(5), xlabel="Wavelength (μm)",ylabel="Scattering cross section",  title="Normalized Scattering Cross Sections",yscale=log10, xscale=log10)
    fill_between!(ax1, [0.76,0.78], [1e-2,1e-2], [9,9], color = (CS[4],0.6))
    fill_between!(ax1, [1.26,1.28], [1e-2,1e-2], [9,9], color = (CS[4],0.6))
    fill_between!(ax1, [1.55,1.75], [1e-2,1e-2], [9,9], color = (CS[5],0.6))
    fill_between!(ax1, [2.04,2.38], [1e-2,1e-2], [9,9], color = (CS[9],0.6))
    text!(ax1,0.76, 1.5, text = L"\text{O_2 A-band}", rotation=π/2, align = (:left,:bottom),fontsize =12)
    text!(ax1,1.25, 1.0, text = L"\text{O_2 1.27\mu\,m band}", rotation=π/2, align = (:left,:bottom),fontsize =12)
    text!(ax1,1.75, 1.0, text = L"\text{CO_2/CH_4 bands}", rotation=π/2, align = (:left,:bottom),fontsize =12)
    text!(ax1,2.38, 0.6, text = L"\text{CO_2/CH_4/N_2O bands}", rotation=π/2, align = (:left,:bottom),fontsize =12)
    #fill_between!(ax1, [0.76,0.78], [1e-2,1e-2], [9,9], color = (CS[4],0.2))
    lines!(ax1, wl/1e3, scattXS[1,1,1,:]./scattXS[1,1,1,n],  label=L"\text{r_g=0.08, \sigma_g=1.3}", color=CarbonI_colors[1], linewidth=3)
    lines!(ax1, wl/1e3, scattXS[1,1,2,:]./scattXS[1,1,2,n],  label=L"\text{r_g=0.08, \sigma_g=1.6}", color=CarbonI_colors[1], linewidth=3, linestyle=:dash)
    lines!(ax1, wl/1e3, scattXS[1,2,1,:]./scattXS[1,2,1,n],  label=L"\text{r_g=0.13, \sigma_g=1.3}", color=CarbonI_colors[7], linewidth=3)
    lines!(ax1, wl/1e3, scattXS[1,2,2,:]./scattXS[1,2,2,n],  label=L"\text{r_g=0.13, \sigma_g=1.6}", color=CarbonI_colors[7], linewidth=3, linestyle=:dash)
    lines!(ax1, wl/1e3, scattXS[1,3,2,:]./scattXS[1,3,2,n],  label=L"\text{r_g=0.20, \sigma_g=1.6}", color=CarbonI_colors[10], linewidth=3, linestyle=:dash)

    #O2:


    CairoMakie.xlims!(ax1, wl[1]/1e3, wl[end]/1e3)
    CairoMakie.ylims!(ax1, 1e-2, 9)
    xticks_values = [0.50, 0.750, 1.0, 1.5000, 2.0000, 2.5000]
    xticks_labels = ["0.5", "0.75", "1.0", "1.5", "2.0", "2.5"]
    ax1.xticks = (xticks_values, xticks_labels)

    ax2 = Axis(f[1,1],xminorgridvisible = true,bottomspinecolor=:gray,leftspinecolor=:gray,backgroundcolor = :transparent, xminorticks = IntervalsBetween(10), spinewidth=2, xlabel="Aerosol radius (μm)",ylabel="Number Density",  title="Microphysical Properties",xscale=log10)
    lines!(ax2, aer_size,pdf(LogNormal(log(0.08), log(1.3)), aer_size),  label=L"\text{r_g=0.08, \sigma_g=1.3}", color=CarbonI_colors[1],linewidth=3)
    lines!(ax2, aer_size,pdf(LogNormal(log(0.08), log(1.6)), aer_size),  label=L"\text{r_g=0.08, \sigma_g=1.6}", color=CarbonI_colors[1], linewidth=3, linestyle=:dash)
    lines!(ax2, aer_size,pdf(LogNormal(log(0.13), log(1.3)), aer_size),  label=L"\text{r_g=0.13, \sigma_g=1.3}", color=CarbonI_colors[7], linewidth=3)
    lines!(ax2, aer_size,pdf(LogNormal(log(0.13), log(1.6)), aer_size),  label=L"\text{r_g=0.13, \sigma_g=1.6}", color=CarbonI_colors[7], linewidth=3, linestyle=:dash)
    lines!(ax2, aer_size,pdf(LogNormal(log(0.2), log(1.6)), aer_size),  label=L"\text{r_g=0.20, \sigma_g=1.6}", color=CarbonI_colors[10], linewidth=3, linestyle=:dash)
    axislegend(ax2, position=:rt, labelsize =12)
    CairoMakie.xlims!(ax2, 0.008, 1.100)
    CairoMakie.ylims!(ax2, 0, 22)
    # Define the desired ticks and labels
    xticks_values = [0.01, 0.10, 1]
    xticks_labels = ["0.01", "0.1", "1"]
    yticks_values = [0.01, 0.1, 1]
    yticks_labels = ["0.01", "0.1", "1.0"]

    # Apply custom tick labels
    ax2.xticks = (xticks_values, xticks_labels)
    ax1.yticks = (yticks_values, yticks_labels)
    hidespines!(ax2, :t, :r) # only top and right

    f
end
f = plotAeroSizes()
save("plots/final/AerosolMicrophysics.pdf", f)

f = with_theme(plotAeroSizes, theme_black())
save("plots/final/AerosolMicrophysics_dark.pdf", f)

function plotAeroSizes_v2()
    f = Figure(resolution=(800,350), title="Aerosol Size Distribution", fontsize=18, backgroundcolor = :transparent, fonts = (; regular = "Helvetica Condensed Light", bold="Helvetica Condensed Bold"))
    ax1 = Axis(f[1,2], yminorgridvisible = true,  yminorticks = IntervalsBetween(5), xlabel="Wavelength (μm)",ylabel="Scattering cross section", backgroundcolor=:transparent,  title="Normalized Scattering Cross Sections",yscale=log10, xscale=log10)
    fill_between!(ax1, [0.76,0.78], [1e-2,1e-2], [9,9], color = (CS[4],0.3))
    fill_between!(ax1, [1.26,1.28], [1e-2,1e-2], [9,9], color = (CS[4],0.3))
    fill_between!(ax1, [1.55,1.75], [1e-2,1e-2], [9,9], color = (CS[5],0.3))
    fill_between!(ax1, [2.04,2.38], [1e-2,1e-2], [9,9], color = (CS[9],0.3))
    text!(ax1,0.76, 1.5, text = L"\text{O_2 A}", rotation=π/2, align = (:left,:bottom),fontsize =15)
    text!(ax1,1.25, 1.0, text = L"\text{O_2 1.27\mu\,m}", rotation=π/2, align = (:left,:bottom),fontsize =15)
    text!(ax1,1.75, 1.0, text = L"\text{CO_2/CH_4}", rotation=π/2, align = (:left,:bottom),fontsize =15)
    text!(ax1,2.38, 0.6, text = L"\text{CO_2/CH_4/N_2O}", rotation=π/2, align = (:left,:bottom),fontsize =15)
    #fill_between!(ax1, [0.76,0.78], [1e-2,1e-2], [9,9], color = (CS[4],0.2))
    lines!(ax1, wl/1e3, scattXS[1,1,1,:]./scattXS[1,1,1,n],  label="Fine Mode", color=CarbonI_colors[4], linewidth=3)
    #lines!(ax1, wl/1e3, scattXS[1,1,2,:]./scattXS[1,1,2,n],  label=L"\text{r_g=0.08, \sigma_g=1.6}", color=CarbonI_colors[1], linewidth=3, linestyle=:dash)
    #lines!(ax1, wl/1e3, scattXS[1,2,1,:]./scattXS[1,2,1,n],  label=L"\text{r_g=0.13, \sigma_g=1.3}", color=CarbonI_colors[7], linewidth=3)
    #lines!(ax1, wl/1e3, scattXS[1,2,2,:]./scattXS[1,2,2,n],  label=L"\text{r_g=0.13, \sigma_g=1.6}", color=CarbonI_colors[7], linewidth=3, linestyle=:dash)
    lines!(ax1, wl/1e3, scattXS[1,3,2,:]./scattXS[1,3,2,n],  label="Coarser Mode", color=CarbonI_colors[5], linewidth=3)

    #O2:


    CairoMakie.xlims!(ax1, wl[1]/1e3, wl[end]/1e3)
    CairoMakie.ylims!(ax1, 1e-2, 9)
    xticks_values = [0.50, 0.750, 1.0, 1.5000, 2.0000, 2.5000]
    xticks_labels = ["0.5", "0.75", "1.0", "1.5", "2.0", "2.5"]
    ax1.xticks = (xticks_values, xticks_labels)

    ax2 = Axis(f[1,1],xminorgridvisible = true,bottomspinecolor=:gray,leftspinecolor=:gray,backgroundcolor=:transparent, xminorticks = IntervalsBetween(10), spinewidth=2, xlabel="Aerosol radius (μm)",ylabel="Number Density",  title="Aerosol Size Distribution",xscale=log10)
    lines!(ax2, aer_size,pdf(LogNormal(log(0.08), log(1.3)), aer_size),  label="Finer Mode", color=CarbonI_colors[4], linewidth=3)
    #lines!(ax2, aer_size,pdf(LogNormal(log(0.13), log(1.3)), aer_size),  label=L"\text{r_g=0.13, \sigma_g=1.3}", color=CarbonI_colors[7], linewidth=3)
    #lines!(ax2, aer_size,pdf(LogNormal(log(0.13), log(1.6)), aer_size),  label=L"\text{r_g=0.13, \sigma_g=1.6}", color=CarbonI_colors[7], linewidth=3, linestyle=:dash)
    lines!(ax2, aer_size,pdf(LogNormal(log(0.2), log(1.6)), aer_size),  label="Coarser Mode", color=CarbonI_colors[5], linewidth=3)
    axislegend(ax2, position=:rt, labelsize =15)
    CairoMakie.xlims!(ax2, 0.008, 1.100)
    CairoMakie.ylims!(ax2, 0, 22)
    # Define the desired ticks and labels
    xticks_values = [0.01, 0.10, 1]
    xticks_labels = ["0.01", "0.1", "1"]
    yticks_values = [0.01, 0.1, 1]
    yticks_labels = ["0.01", "0.1", "1.0"]

    # Apply custom tick labels
    ax2.xticks = (xticks_values, xticks_labels)
    ax1.yticks = (yticks_values, yticks_labels)
    hidespines!(ax2, :t, :r) # only top and right

    f
end
#f = with_theme(plotAeroSizes_v2, theme_ggplot2())
f = plotAeroSizes_v2()
save("plots/final/Box-D3-AerosolMicrophysic_v2.pdf", f)
save("plots/final/Box-D3-AerosolMicrophysic_v2.eps", f)
f = with_theme(plotAeroSizes_v2, theme_black())
save("plots/final/AerosolMicrophysics_dark_v2.pdf", f)

function plotAeroAlbedo()
    f = Figure(resolution=(800,700), title="Aerosol Microphysics", fontsize=18, backgroundcolor = :transparent, fonts = (; regular = "Helvetica Condensed Light", bold="Helvetica Condensed Bold"))
    ax1 = Axis(f[1,2],yminorgridvisible = true,  backgroundcolor = :transparent,yminorticks = IntervalsBetween(5), xlabel="Wavelength (μm)",ylabel="Scattering cross section",  title="Normalized Scattering Cross Sections",yscale=log10, xscale=log10)
    fill_between!(ax1, [0.76,0.78], [1e-2,1e-2], [9,9], color = (CS[4],0.6))
    fill_between!(ax1, [1.26,1.28], [1e-2,1e-2], [9,9], color = (CS[4],0.6))
    fill_between!(ax1, [1.55,1.75], [1e-2,1e-2], [9,9], color = (CS[5],0.6))
    fill_between!(ax1, [2.04,2.38], [1e-2,1e-2], [9,9], color = (CS[9],0.6))
    text!(ax1,0.76, 1.5, text = L"\text{O_2 A-band}", rotation=π/2, align = (:left,:bottom),fontsize =12)
    text!(ax1,1.25, 1.0, text = L"\text{O_2 1.27\mu\,m band}", rotation=π/2, align = (:left,:bottom),fontsize =12)
    text!(ax1,1.75, 1.0, text = L"\text{CO_2/CH_4 bands}", rotation=π/2, align = (:left,:bottom),fontsize =12)
    text!(ax1,2.38, 0.6, text = L"\text{CO_2/CH_4/N_2O bands}", rotation=π/2, align = (:left,:bottom),fontsize =12)
    #fill_between!(ax1, [0.76,0.78], [1e-2,1e-2], [9,9], color = (CS[4],0.2))
    lines!(ax1, wl/1e3, scattXS[1,1,1,:]./scattXS[1,1,1,n],  label=L"\text{r_g=0.08, \sigma_g=1.3}", color=CarbonI_colors[1], linewidth=3)
    lines!(ax1, wl/1e3, scattXS[1,1,2,:]./scattXS[1,1,2,n],  label=L"\text{r_g=0.08, \sigma_g=1.6}", color=CarbonI_colors[1], linewidth=3, linestyle=:dash)
    lines!(ax1, wl/1e3, scattXS[1,2,1,:]./scattXS[1,2,1,n],  label=L"\text{r_g=0.13, \sigma_g=1.3}", color=CarbonI_colors[7], linewidth=3)
    lines!(ax1, wl/1e3, scattXS[1,2,2,:]./scattXS[1,2,2,n],  label=L"\text{r_g=0.13, \sigma_g=1.6}", color=CarbonI_colors[7], linewidth=3, linestyle=:dash)
    lines!(ax1, wl/1e3, scattXS[1,3,2,:]./scattXS[1,3,2,n],  label=L"\text{r_g=0.20, \sigma_g=1.6}", color=CarbonI_colors[10], linewidth=3, linestyle=:dash)

    #O2:


    CairoMakie.xlims!(ax1, wl[1]/1e3, wl[end]/1e3)
    CairoMakie.ylims!(ax1, 1e-2, 9)
    xticks_values = [0.50, 0.750, 1.0, 1.5000, 2.0000, 2.5000]
    xticks_labels = ["0.5", "0.75", "1.0", "1.5", "2.0", "2.5"]
    ax1.xticks = (xticks_values, xticks_labels)

    ax2 = Axis(f[1,1],xminorgridvisible = true,backgroundcolor = :transparent, xminorticks = IntervalsBetween(10), xlabel="Aerosol radius (μm)",ylabel="Number Density",  title="Microphysical Properties",xscale=log10)
    lines!(ax2, aer_size,pdf(LogNormal(log(0.08), log(1.3)), aer_size),  label=L"\text{r_g=0.08, \sigma_g=1.3}", color=CarbonI_colors[1],linewidth=3)
    lines!(ax2, aer_size,pdf(LogNormal(log(0.08), log(1.6)), aer_size),  label=L"\text{r_g=0.08, \sigma_g=1.6}", color=CarbonI_colors[1], linewidth=3, linestyle=:dash)
    lines!(ax2, aer_size,pdf(LogNormal(log(0.13), log(1.3)), aer_size),  label=L"\text{r_g=0.13, \sigma_g=1.3}", color=CarbonI_colors[7], linewidth=3)
    lines!(ax2, aer_size,pdf(LogNormal(log(0.13), log(1.6)), aer_size),  label=L"\text{r_g=0.13, \sigma_g=1.6}", color=CarbonI_colors[7], linewidth=3, linestyle=:dash)
    lines!(ax2, aer_size,pdf(LogNormal(log(0.2), log(1.6)), aer_size),  label=L"\text{r_g=0.20, \sigma_g=1.6}", color=CarbonI_colors[10], linewidth=3, linestyle=:dash)
    axislegend(ax2, position=:rt, labelsize =15)
    CairoMakie.xlims!(ax2, 0.008, 1.100)
    CairoMakie.ylims!(ax2, 0, 22)
    # Define the desired ticks and labels
    xticks_values = [0.01, 0.10, 1]
    xticks_labels = ["0.01", "0.1", "1"]
    yticks_values = [0.01, 0.1, 1]
    yticks_labels = ["0.01", "0.1", "1.0"]

    # Apply custom tick labels
    ax2.xticks = (xticks_values, xticks_labels)
    ax1.yticks = (yticks_values, yticks_labels)

    ax3 = Axis(f[2,2],xminorgridvisible = true,backgroundcolor = :transparent,  xlabel="Wavelength (μm)",ylabel="Albedo",xscale=log10)
    fill_between!(ax3, [0.76,0.78], [0,0], [0.4,0.4], color = (CS[4],0.6))
    fill_between!(ax3, [1.26,1.28], [0,0], [0.4,0.4], color = (CS[4],0.6))
    fill_between!(ax3, [1.55,1.75], [0,0], [0.4,0.4], color = (CS[5],0.6))
    fill_between!(ax3, [2.04,2.38], [0,0], [0.4,0.4], color = (CS[9],0.6))
    lines!(ax3,clima_alb[:,1]/1e3,clima_alb[:,2], color=CS[1], linewidth=3)
    hidexdecorations!(ax1, grid=false)
    CairoMakie.xlims!(ax3, wl[1]/1e3, wl[end]/1e3)
    CairoMakie.ylims!(ax3, 0,0.4)
    xticks_values = [0.50, 0.750, 1.0, 1.5000, 2.0000, 2.5000]
    xticks_labels = ["0.5", "0.75", "1.0", "1.5", "2.0", "2.5"]
    ax3.xticks = (xticks_values, xticks_labels)
    f
end

f = with_theme(plotAeroAlbedo, theme_ggplot2())
save("plots/final/AerosolMicrophysicsAlbedo.pdf", f)
f = with_theme(plotAeroAlbedo, theme_black())
save("plots/final/AerosolMicrophysicsAlbedo_dark.pdf", f)
