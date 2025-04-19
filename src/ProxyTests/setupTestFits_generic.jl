using CarbonI, vSmartMOM
using ImageFiltering, DiffResults, ForwardDiff, InstrumentOperator, Unitful, Interpolations
using NCDatasets, Polynomials, LinearAlgebra, SpecialPolynomials
using CairoMakie, JLD2, Plots, Statistics

"Struct for an atmospheric profile"
struct FitParams{FT}
    "High Resolution Indices (for the forward model)"
    indHR::Array{Int,1}
    "High Resolution Indices (for the forward model)"
    indHR2::Array{Int,1}
    "Low Resolution Indices (for the measurements)"
    indLR::Array{Int,1}
    "Instrument wavelength Arrays"
    wl_ci::Array{FT,1}
    "Convolution Matrix"
    cM::Matrix{FT}
    "Solar Zenith Angle"
    sza::FT
    "Viewing Zenith Angle"
    vza::FT
    "Vertical Column Density (Dry)"
    vcd_dry::Array{FT,1}
    "Cross Sections"
    σ_matrix::Array{FT,3}
    "Gas Profiles"
    gasProfiles::Array{Array{FT,1},1}
    "Atmospheric Profile"
    profile::CarbonI.AtmosphericProfile{FT}
end

# Example Gas Array (we merge CO2_13 and CO2 now) ########################
#=
gas_array = ["co2", "ch4", "n2o", "h2o","hdo"];
#gas_array = ["ch4", "co2", "n2o", "h2o", "co", "hdo", "c2h6"];
#gas_array = ["ch4", "co2", "n2o",  "co",  "c2h6"];
#gas_array = ["ch4", "n2o", "h2o", "co", "hdo", "c2h6"];
setupFile = "/home/cfranken/code/gitHub/CarbonI/src/yaml/carbon-i.yaml"
n_layers = 10
indLR = 8:260
#indLR = 287:410
#indLR = 7:410
cls        = Dict(gas => 1550.0 for gas in gas_array)
cls["h2o"] = 150.0
cls["co2"] = 150.0
cls["hdo"] = 150.0
rel_errors = Dict(gas => 0.15 for gas in gas_array)
rel_errors["h2o"] = 0.15
rel_errors["hdo"] = 0.15
rel_errors["n2o"] = 0.15
pbl_error = 100.0
a = createFitParams(indLR, setupFile, n_layers, gas_array)
x, Sa, h_column = createBayesianConstraints(gas_array, a.gasProfiles,5, a.profile, 150.0, 900.0, cls, rel_errors)
=#
############################################################################

# Need to add options to split profile into strat, free trop and boundary layer
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
    cM, wl_ci = CarbonI.create_carbonI_conv_matrix_cbe(wl)

    # Spectral database is in reverse oder (wavenumber)
    indHR  = findall(wl_ci[indLR[1]]-5 .< reverse(wl) .< wl_ci[indLR[end]]+5)
    indHR2 = findall(wl_ci[indLR[1]]-5 .< wl          .< wl_ci[indLR[end]]+5)

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
    
    # Define prior state vector elements (can be adapted later):
    vmr_co2 = zeros(nL)    .+ parameters.absorption_params.vmr["CO2"] 
    vmr_ch4 = zeros(nL)    .+ parameters.absorption_params.vmr["CH4"]
    vmr_co  = zeros(nL)    .+ parameters.absorption_params.vmr["CO"]
    vmr_n2o = zeros(nL)    .+ parameters.absorption_params.vmr["N2O"]
    vmr_hdo = zeros(nL)    .+ parameters.absorption_params.vmr["HDO"]
    vmr_co2_13 = zeros(nL) .+ parameters.absorption_params.vmr["CO2_13"]
    vmr_c2h6 = zeros(nL)   .+ parameters.absorption_params.vmr["C2H6"]

    # Reduce dimensions, group layers together to get roughly layers of equal pressure difference:
    profile, σ_matrix, indis, gasProfiles = CarbonI.reduce_profile_vcd(n_layers,profile_hr, σ_matrix_hr,[vmr_co2, vmr_ch4, vmr_h2o, vmr_hdo, vmr_n2o, vmr_co,  vmr_c2h6,vmr_co2_13])
    # Merge CO2 and CO2_13:
    σ_matrix[:,:,1] += σ_matrix[:,:,8]
    n_layers = length(indis)

    # Convert cross sections from cm2/molec to m2/mol
    cf = 6.02214076e23/1e4
    return FitParams(indHR, indHR2, collect(indLR), wl_ci, cM, parameters.sza, 0.0, profile.vcd_dry, cf.*σ_matrix[:,:,index_array], gasProfiles[index_array]./cf, profile)
end

function createBayesianConstraints(gas_array, gasProfiles, n_poly,profile, p_strat, p_pbl, correlationLengths, rel_errors, pbl_error)
    n_layers_ = [length(i) for i in gasProfiles];
    n_layers = n_layers_[1] # For now, should enable flexible grid for other gases 
    # Get pressure in hPa from profile
    p = profile.p/1e2;
    # Define a polynomial scaling for the surface polynomial 
    poly = Legendre([1; eps().*ones(n_poly-1)]);
    # Define our state vector (concatenate the gas profiles and polynomial terms):
    x = [reduce(vcat,gasProfiles);poly[:]]
    # Length of the state vector
    n_state = length(x);

    # Create column averaging operators for each gas (need to adapt so that each gas can have a different n_layers):
    h_column = Dict(gas => zeros(n_state) for gas in gas_array)
    for (i,gas) in enumerate(gas_array)
        h_column[gas][(i-1)*n_layers+1:i*n_layers] .= 1.0;
    end
    correlationMatrices = [CarbonI.createCorrelationMatrix(p,p_strat,correlationLengths[gas]) for gas in gas_array];
    errorVectors        = [CarbonI.createErrorVector(p,p_strat,p_pbl, rel_errors[gas], gasProfiles[i]; bl_error=pbl_error) for (i,gas) in enumerate(gas_array)];
    
    # Define prior covariance matrix:
    Sa = zeros(n_state,n_state);
    for (i,gas) in enumerate(gas_array)
        Sa[(i-1)*n_layers+1:i*n_layers,(i-1)*n_layers+1:i*n_layers] .= CarbonI.createPriorCovarianceMatrix(errorVectors[i], correlationMatrices[i]);
    end
    # Put in arbitrarily high numbers for the polynomial term, so these won't be constrained at all:
    counter = 1
    for i=length(gas_array)*n_layers+1:n_state
        if counter < 3
            Sa[i,i] = 1^2;
        else
            Sa[i,i] = 0.001^2;
        end
        counter += 1
    end
    #Sa[n_layers*length(gas_array)+1:n_state,n_layers*length(gas_array)+1:n_state] .= 1e5;

    return x, Sa, h_column
end

# Define the outer function
function define_forward_model(FParams::FitParams{FT}) where {FT}
    (;indHR,indHR2,indLR, wl_ci, cM, sza, vza,σ_matrix ) = FParams
    cM2 = cM[indLR,indHR2]
    wl = wl_ci[indLR]
    σ_matrix_ = σ_matrix[indHR,:,:]
    # x-axis for polynomial [-1,1], enables legendre later:
    @views x_poly = CarbonI.rescale_x(wl)
    # Total sum of τ
    
    dims = size(σ_matrix_)
    # Define the inner function
    function F(x::AbstractArray{FT2}) where{FT2}
        dims = size(σ_matrix)
        vcds = reshape(x[1:(dims[2]*dims[3])],(dims[2],dims[3]) )
        poly = Legendre(x[dims[2]*dims[3]+1:end])
    
        # Air Mass Factor
        AMF = 1/cosd(sza) + 1/cosd(vza);
        # Array for the total sum of τ (adapted to the state vector type)
        ∑τ = zeros(FT2,length(indHR))
        
        for i=1:size(vcds,2)
            for j=1:dims[2]
                ∑τ .+= view(σ_matrix_,:,j,i) * vcds[j,i] 
            end
        end
        # Transmission without Tsolar
        @views T = reverse(exp.(-AMF * ∑τ))
        
        @views T_conv = cM2 * T
        
        #x_poly = CarbonI.rescale_x(wl)
        return T_conv .* poly.(x_poly)
    end

    return F
end


# Define the outer function
function define_inverse_model(F, xa, Sa, Se, n, max_iter,iP ) where {FT}
    result = DiffResults.JacobianResult(zeros(n),xa);
    x_all = zeros(length(xa),max_iter+1)
    K = zeros(n,length(xa));
    # Define the inverse system
    function invert(y::AbstractArray{FT2}) where{FT2}
        # Start at prior (just adjust for first polynomial term)
        x_all[:,1] = xa
        x_all[iP,1] = maximum(y)
        # Start iterations:
        for i=1:max_iter
            ForwardDiff.jacobian!(result, F, x_all[:,i]);
            K .= DiffResults.jacobian(result);
            Fᵢ = DiffResults.value(result);
            iGain = inv(K'inv(Se)K + inv(Sa))K'inv(Se);
            x_all[:,i+1] = xa + iGain * (y - Fᵢ + K *(x_all[:,i]-xa));
        end
        # Posterior covariance matrix:
        Ŝ = inv(K'inv(Se)K + inv(Sa));
        A = (inv(K'inv(Se)K + inv(Sa))K'inv(Se))*K;
        y_mod = F(x_all[:,end]);
        return x_all[:,end], y_mod, Ŝ, A, K
    end
    return invert
end

