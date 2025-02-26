using DelimitedFiles, Plots
# David R Thompson

include("CI_colors.jl")

x = readdlm("data/carboni_cbe_noise.txt")
wl = x[:,1]
h=plot(wl,x[:,2],label="Total Noise",color=CarbonI_colors[1], yaxis=:log, xlabel="Wavelength (nm)", ylabel="Noise, e-", legend=:bottomleft)
h=plot!(wl,x[:,3],label="Shot Noise",color=CarbonI_colors[7])
 #h=plot!(wl,x[:,4],label="Telescope Noise",color=CarbonI_colors[7])
h=plot!(wl,x[:,5],label="Spectrometer Noise",color=CarbonI_colors[10])
h=plot!(wl,x[:,6],label="Quantization Noise",color=CarbonI_colors[13])
h=plot!(wl,x[:,7],label="Dark Noise",color=CarbonI_colors[2])
h=plot!(wl,x[:,8],label="Read Noise",color=CarbonI_colors[5])
plot!(size=(400,400))
savefig("carboni_cbe_noise.svg")
display(h)
 #sleep(5)

x = readdlm("data/carboni_cbe_signal.txt")
wl = x[:,1]
h=plot(wl,x[:,2],label="Scene Signal",color=CarbonI_colors[1], yaxis=:log, xlabel="Wavelength (nm)", ylabel="Signal, e-", legend=:bottomleft)
h=plot!(wl,x[:,3],label="Spectrometer Signal",color=CarbonI_colors[7])
 #h=plot!(wl,x[:,4],label="Telescope Signal",color=CarbonI_colors[7])
h=plot!(wl,x[:,5],label="Dark Signal",color=CarbonI_colors[10])
plot!(size=(400,400))
savefig("carboni_cbe_signal.svg")
display(h)
 #sleep(5)

x = readdlm("data/carboni_cbe_noise.txt")
wl = x[:,1]
h=plot(wl,x[:,2],label="Total Noise",color=CarbonI_colors[1], yaxis=:log, xlabel="Wavelength (nm)", ylabel="Electrons", legend=:bottomleft, linestyle=:dot)
h=plot!(wl,x[:,3],label="Shot Noise",color=CarbonI_colors[7], linestyle=:dot)
 #h=plot!(wl,x[:,4],label="Telescope Noise",color=CarbonI_colors[7])
h=plot!(wl,x[:,5],label="Spectrometer Noise",color=CarbonI_colors[10], linestyle=:dot)
h=plot!(wl,x[:,6],label="Quantization Noise",color=CarbonI_colors[13], linestyle=:dot)
h=plot!(wl,x[:,7],label="Dark Noise",color=CarbonI_colors[2], linestyle=:dot)
h=plot!(wl,x[:,8],label="Read Noise",color=CarbonI_colors[5], linestyle=:dot)
x = readdlm("data/carboni_cbe_signal.txt")
wl = x[:,1]
h=plot!(wl,x[:,2],label="Scene Signal",color=CarbonI_colors[1], yaxis=:log, xlabel="Wavelength (nm)", ylabel="Electrons", legend=:bottomleft)
h=plot!(wl,x[:,3],label="Spectrometer Signal",color=CarbonI_colors[7])
h=plot!(wl,x[:,5],label="Dark Signal",color=CarbonI_colors[10])
plot!(size=(400,400))
savefig("carboni_cbe_signalandnoise.svg")
display(h)
 #sleep(5)


x = readdlm("data/carboni_snr.txt")
wl = x[:,1]
 # estimated
h=plot(wl,x[:,2],color=CarbonI_colors[1], xlabel="Wavelength (nm)", ylabel="Signal to Noise Ratio (SNR)",legend=:false)
 # required
h=plot!(wl,x[:,3],color=CarbonI_colors[7])
annotate!(2130,75,text("Required\nPerformance",8,CarbonI_colors[7]))
annotate!(2280,125,text("Estimated\nPerformance",8,CarbonI_colors[1]))
plot!(size=(400,400))
savefig("carboni_snr.svg")
display(h)
sleep(5)
