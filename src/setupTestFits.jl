
using CarbonI
using ImageFiltering, DiffResults, ForwardDiff, InstrumentOperator, Unitful, Interpolations
using NCDatasets, Polynomials, LinearAlgebra, SpecialPolynomials
using CairoMakie


# Load spectroscopies:
co2, ch4, h2o, hdo, n2o, co, co2_iso2, c2h6 = CarbonI.loadXSModels();

# Read Solar Spectra (not used for the simple tests)
# include("src/readSun.jl")

# Read non-scattering forward model
include("src/forwardModel.jl")

# Define species in the state vector:
hitran_array = (co2, h2o, ch4, co, n2o, hdo, co2_iso2, c2h6);

# Define wavelength grid for forward model:
Δwl = 0.005
wl = 2035:Δwl:2385

parameters = parameters_from_yaml("/home/cfranken/code/gitHub/CarbonI/src/yaml/carbon-i.yaml")

# Generate Atmospheric Profile
profile_hr = CarbonI.generate_atmos_profile(parameters.T, 100*parameters.p, parameters.q)

# Compute cross sections for all layers:
σ_matrix_hr = CarbonI.compute_profile_crossSections(profile_hr, hitran_array , wl);

# Reduce dimensions, group layers together to get roughly layers of equal pressure difference:
n_layers = 10
profile, σ_matrix = CarbonI.reduce_profile(n_layers,profile_hr, σ_matrix_hr)

# Use a flat solar spectrum here (these are small in the 2 micron range anyhow):
solarIrr = ones(length(wl));

# Define an instrument:
FWHM  = 1.0  # 
SSI  = 0.68
kern1 = CarbonI.box_kernel(2*SSI, Δwl)
kern2 = CarbonI.gaussian_kernel(FWHM, Δwl)
kernf = imfilter(kern1, kern2)
lociBox = CarbonI.KernelInstrument(kernf, collect(2040:SSI:2380));
####

# Number of layers:
nL = length(profile.T)
	
# Define prior state vector elements
vmr_co2 = zeros(nL) .+ 400e-6
vmr_ch4 = zeros(nL) .+ 2.0e-6
#vmr_ch4[1:3] .= 1.4e-6
vmr_h2o = zeros(nL) .+0.002 # profile.vcd_h2o ./ profile.vcd_dry
vmr_co  = zeros(nL) .+ 100e-9
vmr_n2o = zeros(nL) .+ 320e-9
#vmr_n2o[1:3] .= 100e-9
vmr_hdo = vmr_h2o * 0.9
vmr_c2h6 = zeros(nL) .+ 1.0e-9
vmrs = [vmr_co2, vmr_h2o, vmr_ch4,vmr_co, vmr_n2o, 0.1*vmr_hdo, 0.1*vmr_co2, 0.1vmr_c2h6];

# Define a polynomial scaling for the surface polynomial
p = Legendre([2.0,0.0001,0.000001,0.0001,0.000001]);

# Define our state vector:
x = [vmr_co2; vmr_h2o; vmr_ch4; vmr_co; vmr_n2o; vmr_hdo;vmr_co2;vmr_c2h6;   p[:] ];

# Get SZA from Parameter file:
sza = parameters.sza

# Get prior covariance matrix:
n_state = length(x);
Sₐ = zeros(n_state,n_state);
rel_error = 0.05;
# vcd_ratio = profile_caltech.vcd_dry ./ mean(profile_caltech.vcd_dry)
	
# Fill the diagonal for the trace gases (hard-coded indices, so we have to be careful here):
for i=1:80
	Sₐ[i,i] = (rel_error*x[i])^2   
end
# higher for H2O
for i=11:19
	Sₐ[i,i] = (0.5*x[i])^2   
end
for i=41:49
	Sₐ[i,i] = (0.5*x[i])^2   
end

# Enlarge prior errors near the surface dramatically:
Sₐ[10,10] = (20*x[10])^2
Sₐ[20,20] = (20*x[20])^2
Sₐ[30,30] = (20*x[30])^2
Sₐ[40,40] = (20*x[40])^2
Sₐ[50,50] = (20*x[50])^2
Sₐ[60,60] = (20*x[60])^2
Sₐ[70,70] = (20*x[70])^2
Sₐ[80,80] = (20*x[80])^2
# Put in arbitrarily high numbers for the polynomial term, so these won't be constrained at all! 
for i=81:n_state
	Sₐ[i,i] = 1e5;
end

# Set the prior state vector:
xa = x;

# Define measurement error covariance matrix:
Se = Diagonal((1e-5* ones(length(lociBox.ν_out))).^2);

#y = CarbonI.conv_spectra(lociBox, wl, spec);


# Start defining column averaging kernels:
h_co2 = zeros(length(x2));
h_co2_ = zeros(length(x2));
h_ch4 = zeros(length(x2));
h_h2o = zeros(length(x2));
h_co  = zeros(length(x2));
h_hdo = zeros(length(x2));
h_n2o = zeros(length(x2));
h_c2h6 = zeros(length(x2));
ratio = profile.vcd_dry/sum(profile.vcd_dry);
h_co2[1:10] .= ratio;
h_h2o[11:20] .= ratio;
h_ch4[21:30] .= ratio;
h_co[31:40] .= ratio;
h_n2o[41:50] .= ratio;
h_hdo[51:60] .= ratio;
h_co2_[61:70] .= ratio;
h_c2h6[71:80] .= ratio;

N = length(lociBox.ν_out)

# Maximum number of iterations:
max_no_of_iter = 6

# Define the state vector for each iteration:
x_all   = zeros((length(x2),max_no_of_iter+1))
F_all   = zeros((N,max_no_of_iter))

# Define measurement here
y = R_conv_carbonI

# Define a polynomial scaling
p = Legendre([maximum(y),0.0001,0.000001,0.0001,0.000001]);
x = [vmr_co2; vmr_h2o; vmr_ch4; vmr_co; vmr_n2o; 0.1*vmr_hdo;0.1*vmr_co2;0.1*vmr_c2h6;   p[:] ];

x_all[:,1]=x


for i=1:max_no_of_iter
    @show i
    #print('Iteration #',i)
	ForwardDiff.jacobian!(result, forward_model_sat_x, x_all[:,i]);
	Kᵢ = DiffResults.jacobian(result);
    Fᵢ = DiffResults.value(result);
    Gain = inv(Kᵢ'inv(Se)Kᵢ + inv(Sₐ))Kᵢ'inv(Se);
    x_all[:,i+1] = xa + Gain * (y - Fᵢ + Kᵢ *(x_all[:,i]-xa))
    #@show h_co2'*x_all[:,end,end]*1e6
    F_all[:,i] = Fᵢ
end

x_final = x_all[:,end]
println("Column averaged CO2: ", h_co2' * x_final * 1e6,"ppm")
println("Column averaged CH4: ", h_ch4' * x_final * 1e9,"ppb")
println("Column averaged N2O: ", h_n2o' * x_final * 1e9,"ppb")
println("Column averaged H2O: ", h_h2o' * x_final * 1e6,"ppm")
println("Column averaged CO: ", h_co' * x_final * 1e9,"ppb")


plot(lociBox.ν_out, F_all[:,end],label="Carbon-I fit")
plot!(lociBox.ν_out, y, label="Carbon-I spectrum (synthetically generated with vSmartMOM)")
