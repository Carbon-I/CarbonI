"Struct for an atmospheric profile"
struct FitParams{FT}
    "High Resolution Indices (for the forward model)"
    indHR::Array{Int,1}
    "High Resolution Indices (for the forward model)"
    indHR2::Array{Int,1}
    "Low Resolution Indices (for the measurements)"
    indLR::Array{Int,1}
    "Instrument wavelength Arrays"
    instrument::Array{FT,1}
    "Convolution Matrix"
    cM::Matrix{FT}
    "Solar Zenith Angle"
    sza::FT
    "Viewing Zenith Angle"
    vza::FT
    "Vertical Column Density (Dry)"
    vcd_dry::Array{FT,1}
    "Cross Sections"
    σ_matrix::Array{Float64,3}
    "Gas Profiles"
    gasProfiles::Array{Array{FT,1},1}
end

gas_array = ["co2", "ch4", "n2o"];


function createFitParams(indLR, setupFile, n_layers, gas_array)
    # Define the dictionary
    name_to_index = Dict(
        "co2" => 1,
        "ch4" => 2,
        "h2o" => 3,
        "hdo" => 4,
        "n2o" => 5,
        "co" => 6,
        "co2_13" => 8,
        "c2h6" => 7
    )

    index_array = [name_to_index[name] for name in gas_array]

    # Load spectroscopies:
    co2, ch4, h2o, hdo, n2o, co, c2h6, co2_iso2 = CarbonI.loadXSModels();

    # Define wavelength grid for forward model:
    Δwl = 0.004
    wl = 2030:Δwl:2390
    # Define an instrument:
    cM, wl_ci = CarbonI.create_carbonI_conv_matrix(wl)

    # Spectral database is in reverse oder (wavenumber)
    indHR = findall(wl_ci[indLR[1]]-5 .< reverse(wl) .< wl_ci[indLR[end]]+5)
    indHR2 = findall(wl_ci[indLR[1]]-5 .< wl .< wl_ci[indLR[end]]+5)

    # Define species in the state vector:
    hitran_array = (co2, ch4, h2o, hdo, n2o, co, co2_iso2, c2h6);

    # Still a fixed yaml input here (could be GeosChem later):
    parameters = parameters_from_yaml(setupFile)

    # Generate Atmospheric Profile
    profile_hr = CarbonI.generate_atmos_profile(parameters.T, 100*parameters.p, parameters.q)

    # Compute cross sections for all layers:
    σ_matrix_hr = CarbonI.compute_profile_crossSections(profile_hr, hitran_array , wl);

    # Trace gas Profiles (make it more abstract later)
    # Number of layers:
    nL = length(profile_hr.T)
    vmr_h2o = parameters.absorption_params.vmr["H2O"]
    # Define prior state vector elements
    vmr_co2 = zeros(nL) .+ parameters.absorption_params.vmr["CO2"] 
    vmr_ch4 = zeros(nL) .+ parameters.absorption_params.vmr["CH4"]
    vmr_co  = zeros(nL) .+ parameters.absorption_params.vmr["CO"]
    vmr_n2o = zeros(nL) .+ parameters.absorption_params.vmr["N2O"]
    vmr_hdo = parameters.absorption_params.vmr["HDO"]
    vmr_co2_13 = zeros(nL) .+ .+ parameters.absorption_params.vmr["CO2_13"]
    vmr_c2h6 = zeros(nL) .+ parameters.absorption_params.vmr["C2H6"]

    # Reduce dimensions, group layers together to get roughly layers of equal pressure difference:
    profile, σ_matrix, indis, gasProfiles = CarbonI.reduce_profile(n_layers,profile_hr, σ_matrix_hr,[vmr_co2, vmr_ch4, vmr_h2o, vmr_hdo, vmr_n2o, vmr_co,  vmr_c2h6,vmr_co2_13])
    # Merge CO2 and CO2_13:
    σ_matrix[:,:,1] += σ_matrix[:,:,8]
    n_layers = length(indis)

    return FitParams(indHR, indHR2, collect(indLR), wl_ci, cM, parameters.sza, 0.0, profile.vcd_dry, σ_matrix[:,:,index_array], gasProfiles[index_array])
end

# Define a polynomial scaling for the surface polynomial
p = Legendre([2.0,0.000001,0.00000001]);

# Define our state vector (concatenate the gas profiles and polynomial terms):
x = [reduce(vcat,gasProfiles);p[:]]

# Get SZA from Parameter file:
sza = parameters.sza

# Get prior covariance matrix:
n_state = length(x);
Sₐ = zeros(n_state,n_state);

rel_error = 0.25;
correlationMatrix = CarbonI.createCorrelationMatrix(profile.p/1e2,0.0,1550)
correlationMatrix_h2o = CarbonI.createCorrelationMatrix(profile.p/1e2,0.0,550)

# Create Diagonal errors for the trace gases:
utls = 250.0
bl   = 850.0
e_ch4  = CarbonI.createErrorVector(profile.p/1e2,utls,bl, rel_error, vmr_ch4)
e_co2  = CarbonI.createErrorVector(profile.p/1e2,utls,bl, rel_error, gasProfiles[1])
#e_co2[end] = vmr_co2[end]*20
e_h2o  = CarbonI.createErrorVector(profile.p/1e2,utls,bl, 2*rel_error, gasProfiles[2])
e_n2o  = CarbonI.createErrorVector(profile.p/1e2,utls,bl, rel_error, gasProfiles[3])
e_co   = CarbonI.createErrorVector(profile.p/1e2,utls,bl, rel_error, vmr_co)
e_hdo  = CarbonI.createErrorVector(profile.p/1e2,utls,bl, 2*rel_error, gasProfiles[4])
e_c2h6 = CarbonI.createErrorVector(profile.p/1e2,utls,bl, rel_error, vmr_c2h6)

# Set center to 0 for Co2 and N2O:
e_co2[2] = 0.001*gasProfiles[1][2]
e_n2o[2] = 0.001*gasProfiles[3][2]

dims = size(σ_matrix)	
# Fill the diagonal for the trace gases (hard-coded indices, so we have to be careful here):
e_coll = (e_co2, e_h2o, e_n2o, e_hdo, e_co2,e_ch4)
iGas = length(e_coll)
index = 1
# Still needs to be smarter for H2O!
for i=1:n_layers:n_layers*iGas
    #@show i,i+n_layers
    if i==n_layers+1
        Sₐ[i:i+n_layers-1,i:i+n_layers-1] .= CarbonI.createPriorCovarianceMatrix(e_coll[index], correlationMatrix_h2o)
    else
        Sₐ[i:i+n_layers-1,i:i+n_layers-1] .= CarbonI.createPriorCovarianceMatrix(e_coll[index], correlationMatrix)
    end
    index += 1
end


# Put in arbitrarily high numbers for the polynomial term, so these won't be constrained at all! 

for i=iGas*n_layers+1:n_state
	Sₐ[i,i] = 1e5;
end

# Set the prior state vector:
xa = x;

# Define measurement error covariance matrix:
errors = 1e10*ones(length(indLR));
errors .= 0.00002
Se = Diagonal(errors.^2);


# Start defining column averaging kernels:
h_co2 = zeros(length(x));
h_co2_ = zeros(length(x));
h_ch4 = zeros(length(x));
h_h2o = zeros(length(x));
h_co  = zeros(length(x));
h_hdo = zeros(length(x));
h_n2o = zeros(length(x));
h_c2h6 = zeros(length(x));
ratio = profile.vcd_dry/sum(profile.vcd_dry);
h_co2[1:n_layers] .= ratio;
h_h2o[1n_layers+1:2n_layers] .= ratio;
h_n2o[2n_layers+1:3n_layers] .= ratio;
h_hdo[3n_layers+1:4n_layers] .= ratio;
h_co2_[4n_layers+1:5n_layers] .= ratio;
#h_c2h6[5n_layers+1:6n_layers] .= ratio;

N = length(indLR)


n_state = length(x);

# Define the state vector for each iteration:
max_no_of_iter = 4
x_all   = zeros((length(x),max_no_of_iter+1))
F_all   = zeros((N,max_no_of_iter))
K = zeros(N,n_state)
A = zeros(n_state,n_state)
result = DiffResults.JacobianResult(zeros(N),x);

x_all[:,1]=x
n2o = []
co2 = []
co2_13 = []
ch4 = []
h2o = []
resi = []
resi2 = []
ys = []
Fs = []
aods = []
albedos = []

sorted_keys = sort(collect(keys(R_conv_carbonI_dict)));

#for i=1:200
for key in sorted_keys
    @show key
    sza = key[1]
    y = R_conv_carbonI_dict[key][indLR];# + 0.0002*randn(N)
    @show size(y)

    # Figure out this index numerically in the future:
    
    push!(ys,y)
    push!(aods,key[3])
    push!(albedos,key[4])

    res = similar(y)
    Ff = similar(y)
    x[iGas*n_layers+1] = maximum(y)
    x_all[:,1]=x
    for i=1:max_no_of_iter
        #@show i
        #print('Iteration #',i)
        ForwardDiff.jacobian!(result, forward_model_sat_x2, x_all[:,i]);
        Kᵢ = DiffResults.jacobian(result);
        Fᵢ = DiffResults.value(result);
        iGain = inv(Kᵢ'inv(Se)Kᵢ + inv(Sₐ))Kᵢ'inv(Se);
        A.=  iGain*Kᵢ
        A.=  iGain*Kᵢ
        K .= Kᵢ
        x_all[:,i+1] = xa + iGain * (y - Fᵢ + Kᵢ *(x_all[:,i]-xa))
        Ff .= Fᵢ
        #@show x_all[:,i+1]*1e6
        #println("Column averaged CO₂: ", (h_co2' * x_all[:,i+1] * 1e6)/400*100)
        #println("Column averaged 13CO₂: ", (h_co2_' * x_all[:,i+1] * 1e6)/400*100)
        if i==max_no_of_iter
            println("Column averaged N₂O: ", (h_n2o' * x_all[:,i+1] * 1e6)/0.32*100)
            println("Column averaged CO₂: ", (h_co2' * x_all[:,i+1] * 1e6)/400*100)
            println("Column averaged 13CO₂: ", (h_co2_' * x_all[:,i+1] * 1e6)/400*100)
            @show x_all[:,i+1]*1e6
            #println("Column averaged CH₄: ", (h_ch4' * x_all[:,i+1] * 1e9)/2000*100)
            #println("Column averaged C2H6: ", (h_c2h6' * x_all[:,i+1] * 1e9)/1.0*100)
            #println("Column averaged H2O: ", (h_h2o' * x_all[:,i+1] * 1e6)/5313.10461*100)
            #println("Column averaged CO: ", (h_co' * x_all[:,i+1] * 1e6)/0.1*100)
            
            #@show x_all[41:50,i+1]
            res .= y-Fᵢ
            
            
            @show maximum(abs.(res))
        end

        #@show 1e6*x_all[1:9,i+1]
        #@show 1e9*x_all[4n_layers+1:5n_layers,i+1]
        #@show x_all[8n_layers+1:n_state,i+1]
        #@show res[ind]' * res[ind]
        #@show h_co2'*x_all[:,end,end]*1e6
        F_all[:,i] = Fᵢ
    end
    push!(Fs, Ff)
    push!(n2o,(h_n2o' * x_all[:,end] * 1e6)/0.32)
    push!(co2,(h_co2' * x_all[:,end] * 1e6)/400)
    push!(co2_13, (h_co2_' * x_all[:,end] * 1e6)/400)
    #push!(ch4,(h_ch4' * x_all[:,end] * 1e9)/2000)
    push!(h2o,(h_h2o' * x_all[:,end] * 1e6)/5313.104611587684)
    #push!(co2_13, (h_co2_' * x_all[:,end] * 1e6)/400)
    push!(resi, sqrt(res' * res))
    push!(resi2, res)
end

n2o_mw1 = deepcopy(n2o)
co2_mw1 = deepcopy(co2)
co2_13_mw1 = deepcopy(co2_13)
h2o_mw1 = deepcopy(h2o)
@save "mw1_fits_all_v2.jld2" n2o_mw1 co2_mw1 co2_13_mw1 h2o_mw1 resi resi2