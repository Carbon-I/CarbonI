using CairoMakie, MakieThemes, DelimitedFiles, InstrumentOperator, Distributions, Interpolations

xx = -0.5:0.002:0.5
FWHM = 0.042
convKernel = create_instrument_kernel(Normal(0, FWHM/2.3), xx)

function convo(x, y, kernel)
    i = sortperm(x)
    λ_grid = x[i]
    F = y[i]
    interp_I = LinearInterpolation(λ_grid, F);
    res = 0.002;
    grid = minimum(λ_grid):res:maximum(λ_grid)
    grid_out = collect(grid[10:10:end-10]);
    FixedKernel = FixedKernelInstrument(kernel, grid_out);
    F_I = interp_I.(grid);
    I_conv = InstrumentOperator.conv_spectra(FixedKernel, grid, F_I)
    return grid_out, I_conv
end

fname0   = "/home/sanghavi/RamanSIFgrid/raylSIF_sza32_alb0p45_psurf1000hpa_rrs_ABO2.dat"
fname0_low   = "/home/sanghavi/RamanSIFgrid/raylSIF_sza32_alb0p15_psurf500hpa_rrs_ABO2.dat"
fname1   = "/home/sanghavi/RamanSIFgrid/rayl_sza32_alb0p45_psurf1000hpa_rrs_ABO2.dat"

fname0_B   = "/home/sanghavi/RamanSIFgrid/raylSIF_sza32_alb0p15_psurf500hpa_rrs_BBO2.dat"
fname1_B   = "/home/sanghavi/RamanSIFgrid/rayl_sza32_alb0p15_psurf500hpa_rrs_BBO2.dat"

specSIF_RRS    = readdlm(fname0)
specSIF_RRS_low    = readdlm(fname0_low)
specnoSIF_RRS  = readdlm(fname1)

specSIF_RRS_B    = readdlm(fname0_B)
specnoSIF_RRS_B    = readdlm(fname1_B)

specnoSIF_RRS  = readdlm(fname1)


wl_abo2 = 1e7./specSIF_RRS[:,1];

wl_bbo2 = 1e7./specSIF_RRS_B[:,1];



preFac = (specSIF_RRS[:,1].^2/1e7)
preFacB = (specSIF_RRS_B[:,1].^2/1e7)

SIF_abo2   = π * preFac .*(specSIF_RRS[:,2]+specSIF_RRS[:,5]-(specnoSIF_RRS[:,2]+specnoSIF_RRS[:,5]))
# Convolved one:
grid_abo2, SIF_abo2_conv   = convo(wl_abo2, specSIF_RRS[:,2]+specSIF_RRS[:,5], convKernel)
grid_abo2, noSIF_abo2_conv = convo(wl_abo2, specnoSIF_RRS[:,2]+specnoSIF_RRS[:,5], convKernel)
SIF_abo2_LR = π * (1e7./grid_abo2).^2/1e7 .* (SIF_abo2_conv-noSIF_abo2_conv)

SIF_bbo2   = π * preFacB .*(specSIF_RRS_B[:,2]+specSIF_RRS_B[:,5]-(specnoSIF_RRS_B[:,2]+specnoSIF_RRS_B[:,5]))
# Convolved one:
grid_bbo2, SIF_bbo2_conv = convo(wl_bbo2, specSIF_RRS_B[:,2]+specSIF_RRS_B[:,5], convKernel)
grid_bbo2, noSIF_bbo2_conv = convo(wl_bbo2, specnoSIF_RRS_B[:,2]+specnoSIF_RRS_B[:,5], convKernel)
SIF_bbo2_LR = π * (1e7./grid_bbo2).^2/1e7 .* (SIF_bbo2_conv-noSIF_bbo2_conv)


Raman_abo2 = preFac .* specSIF_RRS[:,5]
# Convolved one:
grid_abo2, noRaman_abo2_conv = convo(wl_abo2, preFac .* (specSIF_RRS[:,2]), convKernel)
grid_abo2, Raman_abo2_conv = convo(wl_abo2, preFac .* (specSIF_RRS[:,2] + specSIF_RRS[:,5]), convKernel)
Raman_abo2_LR = Raman_abo2_conv - noRaman_abo2_conv

Raman_abo2_low = preFac .* specSIF_RRS_low[:,5]
# Convolved one:
grid_abo2, noRaman_abo2_conv = convo(wl_abo2, preFac .* (specSIF_RRS_low[:,2]), convKernel)
grid_abo2, Raman_abo2_conv = convo(wl_abo2, preFac .* (specSIF_RRS_low[:,2] + specSIF_RRS_low[:,5]), convKernel)
Raman_abo2_low_LR = Raman_abo2_conv - noRaman_abo2_conv

Raman_bbo2 = preFacB .* specSIF_RRS_B[:,5]
# Convolved one:
grid_bbo2, noRaman_bbo2_conv = convo(wl_bbo2, preFacB .* (specSIF_RRS_B[:,2]), convKernel)
grid_bbo2, Raman_bbo2_conv = convo(wl_bbo2, preFacB .* (specSIF_RRS_B[:,2] + specSIF_RRS_B[:,5]), convKernel)
Raman_bbo2_LR = Raman_bbo2_conv - noRaman_bbo2_conv

grid_abo2, Io_conv = convo(wl_abo2, preFac .* specnoSIF_RRS[:,end]/4000, convKernel)
grid_bbo2, Io_conv_B = convo(wl_bbo2, preFacB .* specnoSIF_RRS_B[:,end]/4000, convKernel)

Makie.set_theme!(ggthemr(:fresh))
f = Figure(resolution=(700,430))
ax2 = Axis(f[1,1], halign=:left, title="SIF vs. Raman impact", ylabel="Radiance (W/m²/sr/μm)", xlabel="Wavelength (nm)")
#ax1 = Axis(f[1,2], halign=:right, title="SIF vs. Raman impact", ylabel="Radiance", xlabel="Wavelength (nm)", yaxisposition=:right)

lines!(ax2, wl_abo2, preFac .* specnoSIF_RRS[:,end]/4000,color=:black, alpha=0.25,linewidth=1.0, label="I₀")

lines!(ax2,wl_abo2,SIF_abo2/4,color=:black, linewidth=2, label="TOA SIF")   
#lines!(ax1,1e7./specRRS[:,1],preFac .* specRRS[:,5],color=:red, alpha=0.5,linewidth=1, label="Raman")

wl_all = sort([wl_bbo2; wl_abo2])
#preFacN = (specRRS0_B[:,1].^2/1e7)
lines!(ax2,wl_abo2,Raman_abo2,color=:red, alpha=0.5,linewidth=2, label="Raman; alb=0.45, 1000hPa")
lines!(ax2,wl_abo2,Raman_abo2_low,color=:green, alpha=0.5,linewidth=1, label="Raman; alb=0.15, 500hPa")
lines!(ax2,wl_all,(wl_all/1e3).^-4/12,color=:blue, alpha=0.5,linewidth=1, label="λ⁻⁴")

lines!(ax2,wl_bbo2,SIF_bbo2,color=:black, linewidth=2, label=nothing)
lines!(ax2,wl_bbo2,Raman_bbo2,color=:green, alpha=0.5,linewidth=1, label=nothing)

lines!(ax2, wl_bbo2, preFacB .* specSIF_RRS_B[:,end]/4000,color=:black, alpha=0.25,linewidth=1.0, label=nothing)


axislegend(ax2, position=:rt, framevisible=false, orientation = :horizontal, nbanks = 2)
#axislegend(ax1, position=:lb, framevisible=false)
#ax1.yticklabelsvisible = false

#CairoMakie.xlims!(ax1,746.0,780)
CairoMakie.xlims!(ax2,746.0,780)
CairoMakie.xlims!(ax2,670.0,774)
#CairoMakie.ylims!(ax1,0.0,0.6)
CairoMakie.ylims!(ax2,0.0,0.5)

f


CairoMakie.save("SIF_Raman_HighRes.pdf", f)

#Convolved!
Makie.set_theme!(ggthemr(:fresh))
f = Figure(resolution=(700,430))
ax2 = Axis(f[1,1], halign=:left, title="SIF vs. Raman impact", ylabel="Radiance (W/m²/sr/μm)", xlabel="Wavelength (nm)")
#ax1 = Axis(f[1,2], halign=:right, title="SIF vs. Raman impact", ylabel="Radiance", xlabel="Wavelength (nm)", yaxisposition=:right)

lines!(ax2, grid_abo2, Io_conv ,color=:black, alpha=0.25,linewidth=1.0, label="I₀")
lines!(ax2, grid_bbo2, Io_conv_B,color=:black, alpha=0.25,linewidth=1.0, label=nothing)


lines!(ax2,grid_abo2,SIF_abo2_LR/4,color=:black, linewidth=2, label="TOA SIF")   
#lines!(ax1,1e7./specRRS[:,1],preFac .* specRRS[:,5],color=:red, alpha=0.5,linewidth=1, label="Raman")

wl_all = sort([wl_bbo2; wl_abo2])
#preFacN = (specRRS0_B[:,1].^2/1e7)
lines!(ax2,grid_abo2,Raman_abo2_LR,color=:red, alpha=0.5,linewidth=2, label="Raman; alb=0.45, 1000hPa")
lines!(ax2,grid_abo2,Raman_abo2_low_LR,color=:green, alpha=0.5,linewidth=1, label="Raman; alb=0.15, 500hPa")

lines!(ax2,wl_all,(wl_all/1e3).^-4/12,color=:blue, alpha=0.5,linewidth=1, label="λ⁻⁴")

lines!(ax2,grid_bbo2,SIF_bbo2_LR,color=:black, linewidth=2, label=nothing)
lines!(ax2,grid_bbo2,Raman_bbo2_LR,color=:green, alpha=0.5,linewidth=1, label=nothing)




axislegend(ax2, position=:rt, framevisible=false, orientation = :horizontal, nbanks = 2)
#axislegend(ax1, position=:lb, framevisible=false)
#ax1.yticklabelsvisible = false

#CairoMakie.xlims!(ax1,746.0,780)
CairoMakie.xlims!(ax2,746.0,780)
CairoMakie.xlims!(ax2,670.0,774)
#CairoMakie.ylims!(ax1,0.0,0.6)
CairoMakie.ylims!(ax2,0.0,0.5)

f

CairoMakie.save("SIF_Raman_Conv.pdf", f)