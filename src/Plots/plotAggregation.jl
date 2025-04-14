# Start from BaselineRequirements.jl
using CarbonI
using ImageFiltering, DiffResults, ForwardDiff, InstrumentOperator, Unitful, Interpolations
using NCDatasets, Polynomials, LinearAlgebra, SpecialPolynomials, DelimitedFiles
using CairoMakie
# Load spectroscopies:
co2, ch4, h2o, hdo, n2o, co, co2_iso2, c2h6 = CarbonI.loadXSModels();

#include(joinpath(@__DIR__, "readSun_DC.jl"))
include(joinpath(@__DIR__, "src/readSun.jl"))
function forward_model_x_(𝐱::AbstractArray{FT} ;sun = solarIrr,reflectance=refl, instrument=lociBox, sza=sza, vza=0.0, profile=profile,σ_matrix=σ_matrix, wl=wl) where {FT}
    dims = size(σ_matrix)
	# @show dims
    #xx
    vmrs = reshape(𝐱[1:(dims[2]*dims[3])],(dims[2],dims[3]) )
    poly = Legendre(𝐱[dims[2]*dims[3]+1:end])
    #@show size(vmrs)
    # Air Mass Factor
    AMF = 1/cosd(sza) + 1/cosd(vza)
    #@show sza
    # Total sum of τ
    ∑τ = zeros(FT,size(σ_matrix,1))
	#@show size(vmrs,2)
    for i=1:size(vmrs,2)
        #@show i, vmrs[end,i]
         ∑τ[:] += sum(σ_matrix[:,:,i] .* (vmrs[:,i] .* profile.vcd_dry)', dims=2)
    end
    # Transmission without Tsolar
    T = sun .* reflectance .* reverse(exp.(-AMF * ∑τ))
	#@show T
    T_conv = CarbonI.conv_spectra(instrument, wl, T)
    L = cosd(sza)*T_conv/π;
    # x-axis for polynomial [-1,1], enables legendre later:
    x_poly = CarbonI.rescale_x(instrument.ν_out)
    #@show poly.(x_poly)
   return L .* poly.(x_poly)
end

# Load some profile:
hitran_array = (co2, h2o, ch4, co, n2o, hdo, co2_iso2, c2h6);


# What latitude do we want? Take Caltech as example
myLat = 36.604
myLon = -97.486

myLat = 34.1478
myLon = -118.1445

# Uncomment for high water scenario:
#myLat = 0.0
#myLon = -62
profile_hr = CarbonI.read_atmos_profile_MERRA2(CarbonI.default_merra_file, myLat, myLon, 7);

# Reduce dimensions, group layers together to get roughly layers of equal pressure difference:
n_layers = 10

Δwl = 0.01
wl = 2000:Δwl:2400
σ_matrix_hr = CarbonI.compute_profile_crossSections(profile_hr, hitran_array , wl);
nL = length(profile_hr.T)
	
vmr_co2 = zeros(nL) .+ 407e-6
vmr_ch4 = zeros(nL) .+ 1.8e-6
vmr_ch4[1:3] .= 1.4e-6
vmr_h2o = profile_hr.vcd_h2o ./ profile_hr.vcd_dry
vmr_co  = zeros(nL) .+ 100e-9
vmr_n2o = zeros(nL) .+ 337e-9
vmr_n2o[1:3] .= 100e-9
vmr_hdo = vmr_h2o * 0.9
vmr_c2h6 = zeros(nL) .+ 1.0e-9
vmrs = [vmr_co2, vmr_h2o, vmr_ch4,vmr_co, vmr_n2o, vmr_hdo, vmr_co2, vmr_c2h6];

sol  = CubicSplineInterpolation(range(wlSol[1],wlSol[end], length=length(wlSol)),solar_irr, extrapolation_bc=Interpolations.Flat());
# Reduce to fewer dimensions:
profile, σ_matrix, indis, gasProfiles = CarbonI.reduce_profile(n_layers, profile_hr, σ_matrix_hr,vmrs);



# Define a polynomial scaling
p = Polynomial([0.2,0.0001,0.000001]);

# Define an instrument:
FWHM  = 0.7  # 
SSI  = 0.7
kern1 = CarbonI.box_kernel(2*SSI, Δwl)
kern2 = CarbonI.gaussian_kernel(FWHM, Δwl)
kernf = imfilter(kern1, kern2)
lociBox = CarbonI.KernelInstrument(kernf, collect(2040:SSI:2380));

# Define state vector:
#x = [vmr_co2; vmr_h2o; vmr_ch4; vmr_co; vmr_n2o; vmr_hdo; vmr_co2 ; vmr_c2h6 ;zeros(10) ];
nLeg = 10
xPoly = zeros(nLeg).+eps()
xPoly[1] = 1.0
x = [reduce(vcat,gasProfiles) ; xPoly ];

# Get prior covariance matrix:
n_state = length(x);
Sₐ = zeros(n_state,n_state);
rel_error = 0.01;
# vcd_ratio = profile_caltech.vcd_dry ./ mean(profile_caltech.vcd_dry)
	
# Fill the diagonal for the trace gases:
for i=1:80
	Sₐ[i,i] = (rel_error*x[i])^2   
end
# CO2 at surface, 100% error
Sₐ[10,10] = (20x[10])^2
Sₐ[20,20] = (20x[20])^2
Sₐ[30,30] = (20x[30])^2
Sₐ[40,40] = (22x[40])^2
Sₐ[50,50] = (22x[50])^2
Sₐ[60,60] = (22x[60])^2
Sₐ[70,70] = (122x[70])^2
Sₐ[80,80] = (1022x[80])^2
# Put in arbitrarily high numbers for the polynomial term, so these won't be constrained at all! 
for i=81:n_state
	Sₐ[i,i] = 1e2;
end
ratio = profile.vcd_dry/sum(profile.vcd_dry);

h_co2 = zeros(length(x));
h_co213 = zeros(length(x));
h_ch4 = zeros(length(x));
h_h2o = zeros(length(x));
h_co  = zeros(length(x));
h_hdo = zeros(length(x));
h_n2o = zeros(length(x));
h_c2h6 = zeros(length(x));

h_co2[1:10] .= ratio;
h_h2o[11:20] .= ratio;
h_ch4[21:30] .= ratio;
h_co[31:40] .= ratio;
h_n2o[41:50] .= ratio;
h_hdo[51:60] .= ratio;
h_co213[61:70] .= ratio;
h_c2h6[71:80] .= ratio;

clima_alb = readdlm(CarbonI.albedo_file,',', skipstart=1)
#soil = CubicSplineInterpolation(450:2500,r[:,140], extrapolation_bc=Interpolations.Flat());
soil = CubicSplineInterpolation(300:2400,clima_alb[:,2]/1.16, extrapolation_bc=Interpolations.Flat());
solarIrr = sol(wl);
refl   = soil(wl);

# Define the instrument specs at 400m:
ET  = 44.0u"ms"         # Exposure time
SSI = (2*0.7)u"nm"      # Spectral resolution
Pitch = 18.0u"μm"       # Pixel pitch
FPA_QE = 0.85           # FPA quantum efficiency
Bench_efficiency = 0.65 # Bench efficiency
Fnumber = 2.2           # F-number
readout_noise = 100    # Readout noise
dark_current = 5e3u"1/s" # Dark current
ins = InstrumentOperator.createGratingNoiseModel(ET, Pitch,FPA_QE, Bench_efficiency, Fnumber, SSI, (readout_noise), dark_current);   


albs = 0.03:0.01:0.6
errors = zeros(length(albs), 8)
result = DiffResults.JacobianResult(zeros(length(lociBox.ν_out)),x);
global sza = 30.0
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
co2_error = CubicSplineInterpolation(albs, errors[:,2]/420)


reference_albedos = [0.06, 0.15]

plot(albs, errors[:,5]/330, label="N2O", linewidth=2)

# Let's assume a pixel is 35m across track  and 400m along track (also margin, even though we used buffer in exposure time < 300m along track)

# Define aggregation length scales (in km):
scales = 0.05:0.1:50

# Define cloud fractions:
cloud_fracs = 0.0:0.02:0.99

# Actual pixel area (in global mode)
pixArea = 0.0345 * 0.303 # (per single pixel, not super pixel)
N_global = 12^2/pixArea
N_target = 1^2/(0.035 * 0.035)
# RSS errors
ch4_error_global = sqrt.(n2o_error.(albs).^2 + ch4_error.(albs).^2);
ch4_error_global_2D = ch4_error_global ./ sqrt.((1 .-cloud_fracs).*N_global)' 
ch4_error_target_2D = ch4_error_global ./ sqrt.((1 .-cloud_fracs).*N_target)'

eff_pix_area = zeros(length(scales), length(cloud_fracs))
for (iS,scale) in enumerate(scales)
    for (iC, cloud_frac) in enumerate(cloud_fracs)
        # Compute the effective pixel size:
        eff_pixel_area = scale^2 * (1 - cloud_frac)
        # Store the results in the array:
        eff_pix_area[iS, iC] = eff_pixel_area
    end
end

f = Figure(resolution=(300,300), title="", fontsize=16, backgroundcolor=:transparent)
ax1 = Axis(f[1,1], ylabel="Cloud Fraction", xlabel="Albedo" )
CairoMakie.xlims!(0.03,0.5)
CairoMakie.ylims!(0,0.99)

co = CairoMakie.contourf!(ax1, albs, cloud_fracs, ch4_error_global_2D*1900, levels=[0.2,0.5, 0.75, 1, 2, 3, 5,  8,12],  labels=true, colorrange=(0.01,20), colormap = (:viridis, 0.5), extendhigh = (:orange,0.4), extendlow = (:gray,0.4)); 
CairoMakie.contour!(ax1, albs, cloud_fracs, ch4_error_global_2D*1900, levels=[ 0.75, 1, 2, 3, 5, 8],  labels=true, colorrange=(0.01,20), labelsize = 14,labelcolor=:black, labelfont = :bold); 
save("plots/final/StandardErrors_GlobalMode.pdf", f)
save("plots/final/StandardErrors_GlobalMode.eps", f)
f

f = Figure(resolution=(300,300), title="", fontsize=16, backgroundcolor=:transparent)
ax1 = Axis(f[1,1], ylabel="Cloud Fraction", xlabel="Albedo")#, title="Predicted Standard Error in Target Mode at 1km (in ppb) " )
CairoMakie.xlims!(0.03,0.5)
CairoMakie.ylims!(0,0.99)

co = CairoMakie.contourf!(ax1, albs, cloud_fracs, ch4_error_target_2D*1900, levels=[1, 2, 3, 4,  5, 7.5, 10, 15],  labels=true, colorrange=(1,20), colormap = (:viridis, 0.5), extendhigh = (:orange,0.4), extendlow = (:gray,0.4)); 
CairoMakie.contour!(ax1, albs, cloud_fracs, ch4_error_target_2D*1900, levels=[1, 2, 3, 4, 5,7.5, 10, 15],  labels=true, colorrange=(1,20), labelsize = 14,labelcolor=:black,labelfont = :bold); 
save("plots/final/StandardErrors_TargetMode.pdf", f)
save("plots/final/StandardErrors_TargetMode.eps", f)
f


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
