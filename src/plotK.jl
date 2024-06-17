f = Figure(resolution=(500,400))
ax1 = Axis(f[1,1])
ax2 = Axis(f[2,1] , ylabel="∂F/∂x (%)")
ax3 = Axis(f[3,1],xlabel="Wavelength [nm]")
#ax4 = Axis(f[4,1],xlabel="Wavelength [nm]" )
k2 = 50e2

h2o = lines!(ax2,lociBox.ν_out, k2*K[:,20]*0.00001,label="H₂O, Δ=10ppm", linewidth=1, alpha=0.5)
hdo = lines!(ax2,lociBox.ν_out, k2*K[:,60]*0.00001,label="HDO, Δ=10ppm", linewidth=1, alpha=0.5)
co2 = lines!(ax1,lociBox.ν_out, k2*K[:,10]*2e-6,label="CO₂, Δ=2ppm", linewidth=1)
#co2_ = lines!(ax1,lociBox.ν_out, k2*K[:,70]*vmr_co2[end],label="¹³CO₂")
ch4 = lines!(ax1,lociBox.ν_out, k2*K[:,30]*10e-9,label="CH₄,  Δ=10ppb", linewidth=1)
co  = lines!(ax3,lociBox.ν_out, k2*K[:,40]*50e-9,label="CO,   Δ=50ppb")
n2o = lines!(ax3,lociBox.ν_out, k2*K[:,50]*10e-9,label="N₂O,  Δ=10ppb")
c2h6= lines!(ax3,lociBox.ν_out, k2*K[:,80]*10e-9,label="C₂H₆, Δ=10ppb")
start = 2035
for ax in (ax1,ax2)
    CairoMakie.xlims!(ax,start,2380)
    
    hidexdecorations!(ax, grid=false)
end
axislegend(ax1, position=:cb, framevisible=false)
axislegend(ax2, position=:lb, framevisible=false)
CairoMakie.xlims!(ax3,start,2380)
##CairoMakie.ylims!(ax1,-0.39,0)
#CairoMakie.ylims!(ax2,-0.69,0)
#CairoMakie.ylims!(ax3,-0.015,0)

axislegend(ax3, position=:lb, framevisible=false)
#rowsize!(f.layout,1,Relative(0.7))
#CairoMakie.ylabel!(ax2,"Jacobians (dF/dx)")
rowgap!(f.layout,0)
f
CairoMakie.save("CarbonIJacobian.pdf", f)
CairoMakie.save("CarbonIJacobian.svg", f)

function add_axis_inset(pos=fig[1, 1]; bgcolor=:snow2,
    halign, valign, width=Relative(0.5),height=Relative(0.35),
    alignmode=Mixed(left=5, right=5))
    inset_box = Axis(pos; width, height, halign, valign, alignmode,
        xticklabelsize=12, yticklabelsize=12, backgroundcolor=bgcolor)
    # bring content upfront
    CairoMakietranslate!(inset_box.scene, 0, 0, 10)
    return inset_box
end

1-σ