# Start from BaselineRequirements.jl
using CarbonI
using ImageFiltering, DiffResults, ForwardDiff, InstrumentOperator, Unitful, Interpolations
using NCDatasets, Polynomials, LinearAlgebra, SpecialPolynomials, DelimitedFiles
using CairoMakie
using ImageFiltering
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
FWHM  = 1.7  # 
SSI  = 0.7
kern1 = CarbonI.box_kernel(2*SSI, Δwl)
kern2 = CarbonI.gaussian_kernel(FWHM, Δwl)
kernf = imfilter(kern1, kern2)
lociBox = CarbonI.KernelInstrument(kernf, collect(2040:SSI:2380));

# Define state vector:
#x = [vmr_co2; vmr_h2o; vmr_ch4; vmr_co; vmr_n2o; vmr_hdo; vmr_co2 ; vmr_c2h6 ;zeros(10) ];
nLeg = 18
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
soil = CubicSplineInterpolation(300:2400,clima_alb[:,2]/1.5, extrapolation_bc=Interpolations.Flat());
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


albs = 0.04:0.01:0.6

FWHMs = 0.05:0.03:15.0
errors = zeros(length(FWHMs), 8)
sza = 30.0
for (i,FWHM) in enumerate(FWHMs)
    @show FWHM
    # Define an instrument:
    #FWHM  = 1.7  #
    # 2.7 times oversampling 
    SSI  = FWHM/2.5

    #kern1 = CarbonI.box_kernel(2*SSI, Δwl)
    kern2 = CarbonI.gaussian_kernel(FWHM, Δwl)
    #kernf = imfilter(kern1, kern2)
    ranger = max(2373 - SSI*480, 2035.0) : SSI : 2373
    # Also Adjust polynomial in terms of bandpass, choosing 1 degree per 20nm:
    nLeg = Int(floor((maximum(ranger) - minimum(ranger))/10))+1
    # Put in arbitrarily high numbers for the polynomial term, so these won't be constrained at all! 
    for i=81:n_state
        if i-81 <= nLeg
            Sₐ[i,i] = 1e2^2
        else
            Sₐ[i,i] = 1e-12^2;
        end
    end
    lociBox = CarbonI.KernelInstrument(kern2, collect(ranger));
    ins = InstrumentOperator.createGratingNoiseModel(ET, Pitch,FPA_QE, Bench_efficiency, Fnumber, (2SSI)u"nm", (readout_noise), dark_current);  
    result = DiffResults.JacobianResult(zeros(length(lociBox.ν_out)),x);
    sza = 30
    refl = soil(wl)
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
    @show sqrt(h_ch4' * Ŝ * h_ch4)*1e9
end

 
f = Figure(resolution=(600,300), title="", fontsize=16)
ax1 = Axis(f[1,1], xlabel="FHWM (nm; SSI=FWHM/2.5)", ylabel="Precision (ppm/ppb)",title="Optimal SSI with fixed detector size (bandpass increases with SSI)", yscale=log10 )
#CairoMakie.xlims!(0.04,0.5)
#CairoMakie.ylims!(0,0.94)

CairoMakie.lines!(ax1, FWHMs,errors[:,5], label="N₂O (ppb)",linewidth=2 ); 
CairoMakie.lines!(ax1, FWHMs,errors[:,1], label="CH₄ (ppb)",linewidth=2 ); 
CairoMakie.lines!(ax1, FWHMs,errors[:,2], label="CO₂ (ppm)",linewidth=2 ); 
CairoMakie.lines!(ax1, [2.1,2.1],[1, 20], label="Design Choice",linewidth=4, alpha=0.5, color=:black );
axislegend(ax1; position = :rt, labelsize = 15)
f

ranger = max.(2373 .- FWHMs/3*480, 2035.0)
k = Kernel.gaussian((1.5,))

function plotOptimalRange()
    f = Figure(resolution=(550,400),backgroundcolor = :transparent, title="", fontsize=20,fonts = (; regular = "Helvetica Condensed Light", bold="Helvetica Condensed Bold"))
    ax1 = Axis(f[1,1],xscale=Makie.pseudolog10,xticks = [0.1, 0.4, 0.7, 1, 2, 3,4], xlabel="SSI (nm)", ylabel="Precision (normalized by best)",title="SSI tradeoff for precision with fixed detector size" )
    CairoMakie.xlims!(0.02,1.0)
    CairoMakie.ylims!(0.95,3.5)
    CairoMakie.lines!(ax1, FWHMs/2.5,imfilter(errors[:,5]./minimum(errors[:,5]),k), label="N₂O",linewidth=2 ); 
    CairoMakie.lines!(ax1, FWHMs/2.5,imfilter(errors[:,1]./minimum(errors[:,1]),k), label="CH₄",linewidth=3.5 ); 
    #CairoMakie.lines!(ax1, FWHMs/3,errors[:,2]./minimum(errors[:,2]), label="CO₂",linewidth=2 ); 
    CairoMakie.lines!(ax1, [2.1/3,2.1/3],[0.95, 5], label="Design Choice",linewidth=4, alpha=0.5, color=:gray );
    #ax2 = Axis(f[1,1], ylabel = "Fit Window Start (nm) ")

    #plot2 = lines!(ax2, FWHMs/3, ranger, color = :black)

    #ax2.yaxisposition = :right
    #ax2.yticklabelalign = (:left, :center)
    #ax2.xticklabelsvisible = false
    ##ax2.xticklabelsvisible = false
    #ax2.xlabelvisible = false

    #linkxaxes!(ax1,ax2)
    CairoMakie.xlims!(0.02,4.0)


    axislegend(ax1; position = :rt, labelsize = 18)
    
    f
end
f = with_theme(plotOptimalRange, theme_ggplot2());
save("plots/final/OptimizedRange.pdf", f)
save("plots/final/Box-D6-OptimizedRange.eps", f)
 f = with_theme(plotOptimalRange, theme_black())
save("plots/final/OptimizedRange_dark.pdf", f)

