using CairoMakie. ColorSchemes

vod = zeros(size(σ_matrix,1), size(σ_matrix,3))
for i in 1:8
    for j=1:10
        vod[:,i] += σ_matrix[:,j,i] * gasProfiles[i][j] * profile.vcd_dry[j]
    end
end
vod[vod.<1e-10] .= 1e-10

# Define an instrument:
FWHM  = 0.6  # 
SSI  = 0.7
kern1 = CarbonI.box_kernel(2*SSI, Δwl)
kern2 = CarbonI.gaussian_kernel(FWHM, Δwl)
kernf = imfilter(kern1, kern2)
CI_Box = CarbonI.KernelInstrument(kernf, collect(2040:SSI:2380));

# Define an instrument:
FWHM  = 8  # 
SSI  = 7
kern1 = CarbonI.box_kernel(FWHM, Δwl)
#kern2 = CarbonI.gaussian_kernel(FWHM, Δwl)
#kernf = imfilter(kern1, kern2)
emit_Box = CarbonI.KernelInstrument(kernf, collect(2040:SSI:2380));

# Define an instrument:
FWHM  = 0.24  # 
SSI  = 0.1
kern2 = CarbonI.gaussian_kernel(FWHM, Δwl)
tropoBox = CarbonI.KernelInstrument(kernf, collect(2040:SSI:2380));

#CarbonI.conv_spectra(instrument, wl, T)

labels = ["CO₂", "H₂O", "CH₄","CO", "N₂O", "HDO", "¹³CO₂", "C₂H₆"]

f = Figure(resolution=(700,400), title="Vertical Optical Depth of Trace Gases", fontsize=16)
ax1 = Axis(f[1,1], ylabel="Optical Depth",  xlabel="Wavelength (nm)", yscale=log10) 
for i in 1:8
    lines!(ax1, CI_Box.ν_out, color=ColorSchemes.tab10[i], -log.(CarbonI.conv_spectra(CI_Box, wl, exp.(-reverse(vod[:,i])))),  linewidth=3, alpha=0.6, label=labels[i])
end
CairoMakie.xlims!(2040,2380)
CairoMakie.ylims!(1e-5,1)
#color=ColorSchemes.Set3_8.colors,
#axislegend(ax1, position=:rt)
f[1, 2] = Legend(f, ax1, "Trace Gases")
f
save("plots/Carbon-I_VODs.pdf", f)

f = Figure(resolution=(700,400), title="Vertical Optical Depth of Trace Gases", fontsize=16)
ax1 = Axis(f[1,1], ylabel="Optical Depth",  xlabel="Wavelength (nm)", yscale=log10) 
for i in 1:8
    lines!(ax1, emit_Box.ν_out, color=ColorSchemes.tab10[i], -log.(CarbonI.conv_spectra(emit_Box, wl, exp.(-reverse(vod[:,i])))),  linewidth=3, alpha=0.6, label=labels[i])
end
CairoMakie.xlims!(2040,2380)
CairoMakie.ylims!(1e-5,1)
#color=ColorSchemes.Set3_8.colors,
#axislegend(ax1, position=:rt)
f[1, 2] = Legend(f, ax1, "Trace Gases")
f
save("plots/EMIT_VODs.pdf", f)
