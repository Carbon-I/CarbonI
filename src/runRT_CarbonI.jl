using vSmartMOM, InstrumentOperator, CarbonI, ImageFiltering, Interpolations, Statistics, CUDA
# device!(1)
Δwl = 0.005
wl = 2035:Δwl:2385
FWHM  = 1.0  # 
SSI  = 0.68
kern1 = CarbonI.box_kernel(2*SSI, Δwl)
kern2 = CarbonI.gaussian_kernel(FWHM, Δwl)
kernf = imfilter(kern1, kern2)
wl_loci = collect(2040:SSI:2380)
lociBox = CarbonI.KernelInstrument(kernf, wl_loci);

wl_emit = collect(2040:7.5:2380)
emitBox = CarbonI.KernelInstrument(CarbonI.box_kernel(8.5, Δwl), wl_emit);

parameters = parameters_from_yaml("/home/cfranken/code/gitHub/CarbonI/src/yaml/carbon-i.yaml")
model = model_from_parameters(parameters);
a = rt_run(model)

aod_band = sum(model.τ_aer[1]); # Total AOD

ν = parameters.spec_bands[1];
#Δν = mean(diff(ν));
#grid = ν[1]:Δν:ν[end]+10eps()
grid = LinRange(ν[1], ν[end], length(ν))

# Get reflected radiance:
R = a[1][1,1,:];

# Interpolate to the grid:
rad_inter = CubicSplineInterpolation(grid, R);
R_conv_carbonI = CarbonI.conv_spectra(lociBox, wl, rad_inter(1e7./wl))
R_conv_emit = CarbonI.conv_spectra(emitBox, wl, rad_inter(1e7./wl))


# Run medium AOD case:
parameters = parameters_from_yaml("src/yaml/carbon-i-mediumAOD.yaml")
model = model_from_parameters(parameters);
a2 = rt_run(model)
rad_inter2 = CubicSplineInterpolation(grid, a2[1][1,1,:]);
R_conv2 = CarbonI.conv_spectra(lociBox, wl, rad_inter2(1e7./wl))


# Run N2O case:
parameters = parameters_from_yaml("src/yaml/carbon-i-dN2O.yaml")
model = model_from_parameters(parameters);
a3 = rt_run(model)
rad_inter3 = CubicSplineInterpolation(grid, a3[1][1,1,:]);
R_conv3 = CarbonI.conv_spectra(lociBox, wl, rad_inter3(1e7./wl))

#max_no_of_iter = 3

#N = length(lociBox.ν_out)
#x[81] = 0.0005
#x[82:end] .= 0.00001
#x_all   = zeros((length(x),max_no_of_iter+1))
#F_all   = zeros((N,max_no_of_iter))
#x_all[:,1]=x

#y = R_conv3
#for i=1:max_no_of_iter
#    @show i
    #print('Iteration #',i)
#	ForwardDiff.jacobian!(result, forward_model_x, x_all[:,i]);#
#	Kᵢ = DiffResults.jacobian(result);
#    Fᵢ = DiffResults.value(result);
#    Gain = inv(Kᵢ'inv(Se)Kᵢ + inv(Sₐ))Kᵢ'inv(Se);
#    x_all[:,i+1] = xa + Gain * (y - Fᵢ + Kᵢ *(x_all[:,i]-xa))
    #@show h_co2'*x_all[:,end,end]*1e6
#    F_all[:,i] = Fᵢ
#end










