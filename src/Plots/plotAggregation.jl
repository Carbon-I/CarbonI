# Start from BaselineRequirements.jl
# Define the instrument specs at 400m:
ET  = 40.0u"ms"         # Exposure time
SSI = (2*0.7)u"nm"      # Spectral resolution
Pitch = 18.0u"μm"       # Pixel pitch
FPA_QE = 0.85           # FPA quantum efficiency
Bench_efficiency = 0.65 # Bench efficiency
Fnumber = 2.2           # F-number
readout_noise = 100    # Readout noise
dark_current = 10e3u"1/s" # Dark current
ins = InstrumentOperator.createGratingNoiseModel(ET, Pitch,FPA_QE, Bench_efficiency, Fnumber, SSI, (readout_noise), dark_current);   


albs = 0.04:0.01:0.6
errors = zeros(length(albs), 8)

for (i,alb) in enumerate(albs)
    sza = 30
    refl = alb.+0.0*soil(wl)
    ForwardDiff.jacobian!(result, forward_model_x_, x);
    K = DiffResults.jacobian(result);
    F = DiffResults.value(result);
    nesr = InstrumentOperator.noise_equivalent_radiance(ins, (lociBox.ν_out)u"nm", (F)u"mW/m^2/nm/sr");
    nesr_ = nesr./1u"mW/m^2/nm/sr"
    Se = Diagonal(nesr_.^2);
    Ŝ = inv(K'inv(Se)K + inv(Sₐ));
    errors[i,1] = sqrt(h_ch4' * Ŝ * h_ch4)*1e9
    errors[i,2] = sqrt(h_co2' * Ŝ * h_co2)*1e6
    errors[i,3] = sqrt(h_h2o' * Ŝ * h_h2o)*1e6
    errors[i,4] = sqrt(h_hdo' * Ŝ * h_hdo)*1e6
    errors[i,5] = sqrt(h_n2o' * Ŝ * h_n2o)*1e9
    errors[i,6] = sqrt(h_co'  * Ŝ * h_co)*1e9
    errors[i,7] = sqrt(h_co213'  * Ŝ * h_co213)*1e6
    errors[i,8] = sqrt(h_c2h6' * Ŝ * h_c2h6)*1e9
end

 
n2o_error = CubicSplineInterpolation(albs, errors[:,5]/330)
ch4_error = CubicSplineInterpolation(albs, errors[:,1]/1900)
co2_error = CubicSplineInterpolation(albs, errors[:,2]/410)


reference_albedos = [0.06, 0.15]

plot(albs, errors[:,5]/330, label="N2O", linewidth=2)

# Let's assume a pixel is 35m across track (i.e. has margin) and 400m along track (also margin, even though we used buffer in exposure time < 300m along track)

# Define aggregation length scales (in km):
scales = 0.05:0.1:50

# Define cloud fractions:
cloud_fracs = 0.0:0.02:0.95

# Actual pixel area (in global mode)
pixArea = 0.035 * 0.4

eff_pix_area = zeros(length(scales), length(cloud_fracs))
for (iS,scale) in enumerate(scales)
    for (iC, cloud_frac) in enumerate(cloud_fracs)
        # Compute the effective pixel size:
        eff_pixel_area = scale^2 * (1 - cloud_frac)
        # Store the results in the array:
        eff_pix_area[iS, iC] = eff_pixel_area
    end
end

custom_label = x -> string(round(x; digits=3), "%")

f = Figure(resolution=(600,240), title="Aggregation length scales when using the N₂O proxy for high accuracy (precision in %)", fontsize=16)
ax1 = Axis(f[1,1], ylabel="Cloud Fraction",  title="Albedo=0.05", )
CairoMakie.xlims!(0,50)
CairoMakie.ylims!(0,0.94)
ax2 = Axis(f[1,2], ylabel="Cloud Fraction", xlabel="Aggregation Scale (km)", title="Albedo=0.15")
CairoMakie.xlims!(0,25)
CairoMakie.ylims!(0,0.94)
ax3 = Axis(f[1,3], ylabel="Cloud Fraction",  title="Albedo=0.3")    # Adjust title font size
    # Adjust x and y axis label font size)
CairoMakie.xlims!(0,15)
CairoMakie.ylims!(0,0.94)
iS = findall(scales .< 50)
co = CairoMakie.contourf!(ax1, scales[iS], cloud_fracs,1.0./sqrt.(eff_pix_area[iS,:]./pixArea)*n2o_error(0.05)*100, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), colormap = (:viridis, 0.5), extendhigh = (:orange,0.4), extendlow = (:gray,0.4)); 
CairoMakie.contour!(ax1, scales[iS ], cloud_fracs,1.0./sqrt.(eff_pix_area[iS,:]./pixArea)*n2o_error(0.05)*100, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), labelsize = 14,labelfont = :bold); 
iS = findall(scales .< 25)
CairoMakie.contourf!(ax2, scales[iS ], cloud_fracs,1.0./sqrt.(eff_pix_area[iS,:]./pixArea)*n2o_error(0.15)*100, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), colormap = (:viridis, 0.5), extendhigh = (:orange,0.4), extendlow = (:gray,0.4)); 
CairoMakie.contour!(ax2, scales[iS ], cloud_fracs,1.0./sqrt.(eff_pix_area[iS,:]./pixArea)*n2o_error(0.15)*100, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), labelsize = 14,labelfont = :bold); 
iS = findall(scales .< 15)
CairoMakie.contourf!(ax3, scales[iS ], cloud_fracs,1.0./sqrt.(eff_pix_area[iS,:]./pixArea)*n2o_error(0.3)*100, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), colormap = (:viridis, 0.5), extendhigh = (:orange,0.4), extendlow = (:gray,0.4)); 
CairoMakie.contour!(ax3, scales[iS ], cloud_fracs,1.0./sqrt.(eff_pix_area[iS,:]./pixArea)*n2o_error(0.3)*100, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), labelsize = 14,labelfont = :bold); 
#CairoMakie.Colorbar(f[1,4],co)
hideydecorations!(ax2, grid=false)
hideydecorations!(ax3, grid=false)

f
save("plots/aggregation_length_scales_GlobalMode.pdf", f)

f = Figure(resolution=(400,600), title="Aggregation length scales when using the N₂O proxy for high accuracy (precision in %)", fontsize=16)
ax1 = Axis(f[1,1], ylabel="Cloud Fraction",  title="Albedo=0.05", )
CairoMakie.xlims!(0,20)
CairoMakie.ylims!(0,0.94)
ax2 = Axis(f[2,1], ylabel="Cloud Fraction", xlabel="Aggregation Scale (km)", title="Albedo=0.15")
CairoMakie.xlims!(0,20)
CairoMakie.ylims!(0,0.94)

iS = findall(scales .< 20)
co = CairoMakie.contourf!(ax1, scales[iS], cloud_fracs,1.0./sqrt.(eff_pix_area[iS,:]./pixArea)*n2o_error(0.05)*100, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), colormap = (:viridis, 0.5), extendhigh = (:orange,0.4), extendlow = (:gray,0.4)); 
CairoMakie.contour!(ax1, scales[iS ], cloud_fracs,1.0./sqrt.(eff_pix_area[iS,:]./pixArea)*n2o_error(0.05)*100, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), labelsize = 14,labelfont = :bold); 
#S = findall(scales .< 25)
CairoMakie.contourf!(ax2, scales[iS ], cloud_fracs,1.0./sqrt.(eff_pix_area[iS,:]./pixArea)*n2o_error(0.15)*100, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), colormap = (:viridis, 0.5), extendhigh = (:orange,0.4), extendlow = (:gray,0.4)); 
CairoMakie.contour!(ax2, scales[iS ], cloud_fracs,1.0./sqrt.(eff_pix_area[iS,:]./pixArea)*n2o_error(0.15)*100, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), labelsize = 14,labelfont = :bold); 

#CairoMakie.Colorbar(f[1,4],co)
#hideydecorations!(ax2, grid=false)
hidexdecorations!(ax1, grid=false)
#hideydecorations!(ax2, grid=false)

f
save("plots/aggregation_length_scales_GlobalMode_2P.pdf", f)



f = Figure(resolution=(400,300), title="Aggregation length scales when using the N₂O proxy for high accuracy (precision in %)")
ax1 = Axis(f[1,1], ylabel="Cloud Fraction", xlabel="Aggregation Scale (km)", title="Albedo=0.1")
CairoMakie.xlims!(0,50)
CairoMakie.ylims!(0,0.94)
#ax2 = Axis(f[1,2], ylabel="Cloud Fraction", xlabel="Aggregation Scale (km)", title="Albedo=0.15")
#CairoMakie.xlims!(0,25)
#CairoMakie.ylims!(0,0.94)
ax3 = Axis(f[1,2], ylabel="Cloud Fraction", xlabel="Aggregation Scale (km)", title="Albedo=0.3")
CairoMakie.xlims!(0,10)
CairoMakie.ylims!(0,0.94)
CairoMakie.contourf!(ax1, scales, cloud_fracs,1.0./sqrt.(eff_pix_area./pixArea)*errors[7,5]/3.30, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), colormap = (:viridis, 0.5), extendhigh = (:orange,0.4), extendlow = (:gray,0.4)); 
CairoMakie.contour!(ax1, scales, cloud_fracs,1.0./sqrt.(eff_pix_area./pixArea)*errors[7,5]/3.30, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), labelsize = 13,labelfont = :bold); 
#CairoMakie.contourf!(ax2, scales, cloud_fracs,1.0./sqrt.(eff_pix_area./pixArea)*4.1, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), colormap = (:viridis, 0.5), extendhigh = (:orange,0.4), extendlow = (:gray,0.4)); 
#CairoMakie.contour!(ax2, scales, cloud_fracs,1.0./sqrt.(eff_pix_area./pixArea)*4.1, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), labelsize = 13,labelfont = :bold); 
CairoMakie.contourf!(ax3, scales, cloud_fracs,1.0./sqrt.(eff_pix_area./pixArea)*errors[27,5]/3.30, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), colormap = (:viridis, 0.5), extendhigh = (:orange,0.4), extendlow = (:gray,0.4)); 
CairoMakie.contour!(ax3, scales, cloud_fracs,1.0./sqrt.(eff_pix_area./pixArea)*errors[27,5]/3.30, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), labelsize = 13,labelfont = :bold)
#CairoMakie.Colorbar(f[1,4],co)
#ideydecorations!(ax2, grid=false)
hideydecorations!(ax3, grid=false)

f
save("plots/aggregation_length_scales_GlobalMode.pdf", f)

# Actual pixel area (in global mode)
pixAreaTarget = 0.035 * 0.035
f = Figure(resolution=(600,240), title="Aggregation length scales when using the N₂O proxy for high accuracy (precision in %)", fontsize=16)
ax1 = Axis(f[1,1], ylabel="Cloud Fraction",  title="Albedo=0.05")
CairoMakie.xlims!(0,14)
CairoMakie.ylims!(0,0.94)
ax2 = Axis(f[1,2], ylabel="Cloud Fraction", xlabel="Aggregation Scale (km)", title="Albedo=0.15")
CairoMakie.xlims!(0,7.5)
CairoMakie.ylims!(0,0.94)
ax3 = Axis(f[1,3], ylabel="Cloud Fraction",  title="Albedo=0.3")
CairoMakie.xlims!(0,5)
CairoMakie.ylims!(0,0.94)
iS = findall(scales .< 14)
co = CairoMakie.contourf!(ax1, scales[iS], cloud_fracs,1.0./sqrt.(eff_pix_area[iS,:]./pixAreaTarget)*n2o_error(0.05)*100, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), colormap = (:viridis, 0.5), extendhigh = (:orange,0.4), extendlow = (:gray,0.4)); 
co1 = CairoMakie.contour!(ax1, scales[iS], cloud_fracs,1.0./sqrt.(eff_pix_area[iS,:]./pixAreaTarget)*n2o_error(0.05)*100, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, label=custom_label, colorrange=(0.05,5), labelsize = 14,labelfont = :bold); 
iS = findall(scales .< 8)
CairoMakie.contourf!(ax2, scales[iS], cloud_fracs,1.0./sqrt.(eff_pix_area[iS,:]./pixAreaTarget)*n2o_error(0.15)*100, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), colormap = (:viridis, 0.5), extendhigh = (:orange,0.4), extendlow = (:gray,0.4)); 
CairoMakie.contour!(ax2, scales[iS], cloud_fracs,1.0./sqrt.(eff_pix_area[iS,:]./pixAreaTarget)*n2o_error(0.15)*100, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), labelsize = 14,labelfont = :bold); 
iS = findall(scales .< 5)
CairoMakie.contourf!(ax3, scales[iS], cloud_fracs,1.0./sqrt.(eff_pix_area[iS,:]./pixAreaTarget)*n2o_error(0.3)*100, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), colormap = (:viridis, 0.5), extendhigh = (:orange,0.4), extendlow = (:gray,0.4)); 
CairoMakie.contour!(ax3, scales[iS], cloud_fracs,1.0./sqrt.(eff_pix_area[iS,:]./pixAreaTarget)*n2o_error(0.3)*100, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), labelsize = 14,labelfont = :bold); 
#CairoMakie.Colorbar(f[1,4],co)
hideydecorations!(ax2, grid=false)
hideydecorations!(ax3, grid=false)

f
save("plots/aggregation_length_scales_TargetMode.pdf", f)

pixAreaTarget = 0.035 * 0.035
f = Figure(resolution=(400,600), title="Aggregation length scales when using the N₂O proxy for high accuracy (precision in %)", fontsize=16)
ax1 = Axis(f[1,1], ylabel="Cloud Fraction",  title="Albedo=0.05")
CairoMakie.xlims!(0,20)
CairoMakie.ylims!(0,0.94)
ax2 = Axis(f[2,1], ylabel="Cloud Fraction", xlabel="Aggregation Scale (km)", title="Albedo=0.15")
CairoMakie.xlims!(0,20)
CairoMakie.ylims!(0,0.94)

iS = findall(scales .< 20)
co = CairoMakie.contourf!(ax1, scales[iS], cloud_fracs,1.0./sqrt.(eff_pix_area[iS,:]./pixAreaTarget)*n2o_error(0.05)*100, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), colormap = (:viridis, 0.5), extendhigh = (:orange,0.4), extendlow = (:gray,0.4)); 
co1 = CairoMakie.contour!(ax1, scales[iS], cloud_fracs,1.0./sqrt.(eff_pix_area[iS,:]./pixAreaTarget)*n2o_error(0.05)*100, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true,  colorrange=(0.05,5), labelsize = 14,labelfont = :bold); 

CairoMakie.contourf!(ax2, scales[iS], cloud_fracs,1.0./sqrt.(eff_pix_area[iS,:]./pixAreaTarget)*n2o_error(0.15)*100, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), colormap = (:viridis, 0.5), extendhigh = (:orange,0.4), extendlow = (:gray,0.4)); 
CairoMakie.contour!(ax2, scales[iS], cloud_fracs,1.0./sqrt.(eff_pix_area[iS,:]./pixAreaTarget)*n2o_error(0.15)*100, levels=[0.05, 0.1,0.15, 0.25, 1],  labels=true, colorrange=(0.05,5), labelsize = 14,labelfont = :bold); 

#CairoMakie.Colorbar(f[1,4],co)
hidexdecorations!(ax1, grid=false)

f
save("plots/aggregation_length_scales_TargetMode_2P.pdf", f)
