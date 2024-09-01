using NCDatasets, Plots, Polynomials

# Compare paper https://acp.copernicus.org/articles/19/9371/2019/#&gid=1&pid=1

oco2 = Dataset("data/oco2_LtCO2_180825_B11100Ar_230602034736s.nc4");
ind = findall(52.5 .< lat .< 53 .&& 38.5 .< lon .< 40.25)

lat = oco2["latitude"][:]
lon = oco2["longitude"][:]
sza = oco2["solar_zenith_angle"][:] 
xco2_s = oco2.group["Preprocessors"]["xco2_strong_idp"][:]
xco2_w = oco2.group["Preprocessors"]["xco2_weak_idp"][:]
xco2_raw = oco2.group["Retrieval"]["xco2_raw"][:]
xco2 = oco2["xco2"][:]
albedo_s = oco2.group["Retrieval"]["albedo_sco2"][:]
albedo_w = oco2.group["Retrieval"]["albedo_wco2"][:]

p = fit(convert.(Float32,albedo_s[ind]), convert.(Float32,xco2_s[ind]), 2);

xco2_s_corr = xco2_s[ind].- p.(convert.(Float32,albedo_s[ind]))

# Plot Albedo:
scatter(albedo_s[ind], xco2_s[ind], label="Strong CO₂ Non-Scattering")
scatter!(albedo_s[ind], p.(convert.(Float32,albedo_s[ind])), label="Polynomial fit")
xlabel!("Albedo")
ylabel!("XCO₂")

# Plot XCO2:
scatter(xco2_raw[ind].-405, xco2_s_corr, label="Raw FP")
scatter!(xco2[ind].-405, xco2_s_corr, label="Bias corrected FP")
xlabel!("Full Physics XCO₂ - 405 ppm")
ylabel!("Albedo-corrected non-scattering XCO₂")


# Map them:
scatter(lon[ind], lat[ind], zcolor=xco2[ind], label="XCO2 FP bias-corrected")
scatter(lon[ind], lat[ind], zcolor=xco2_raw[ind], label="XCO2 FP raw")
scatter(lon[ind], lat[ind], zcolor=xco2_s[ind], label="XCO2 Non-scattering")
scatter!(lon[ind], lat[ind], zcolor=xco2_s_corr, label="XCO2 Non-scattering bias-corrected")

# From silumations:
using JLD2
@load "simulated_rads_all.jld2" R_conv_carbonI_dict
@load "mw2_fits_all.jld2"
@load "mw1_fits_all.jld2"

sorted_keys = sort(collect(keys(R_conv_carbonI_dict)));

paods = [a[2] for a in sorted_keys];
szas = [a[1] for a in sorted_keys];

ind_lowlowAOD = findall(paods .== 850 .&& szas .==40 .&& 0.007 .< aods .< 0.01 .&& albedos .> 0.035 )
ind_lowAOD = findall(paods .== 850 .&& szas .==40 .&& 0.04 .< aods .< 0.06 .&& albedos .> 0.035 )
ind_highAOD = findall(paods .== 850 .&& szas .==40 .&& 0.16 .< aods .< 0.2 .&& albedos .> 0.035) 

scatter(albedos[ind_lowlowAOD], n2o_mw2[ind_lowlowAOD], label="XN₂O; AOD = 0.01")
scatter!(albedos[ind_lowAOD], n2o_mw2[ind_lowAOD], label="XN₂O; AOD = 0.05")
scatter!(albedos[ind_highAOD], n2o_mw2[ind_highAOD], label="XN₂O; AOD = 0.18")
xlabel!("Albedo")
ylabel!("Retrieved/True Ratio")

scatter(albedos[ind_lowlowAOD], ch4_mw2[ind_lowlowAOD], label="XCH₄; AOD = 0.01")
scatter!(albedos[ind_lowAOD], ch4_mw2[ind_lowAOD], label="XCH₄; AOD = 0.05")
scatter!(albedos[ind_highAOD], ch4_mw2[ind_highAOD], label="XCH₄; AOD = 0.18")
xlabel!("Albedo")
ylabel!("Retrieved/True Ratio")

scatter(albedos[ind_lowlowAOD], ch4_mw2[ind_lowlowAOD]./n2o_mw2[ind_lowlowAOD], label="XCH₄/XN₂O; AOD = 0.01", ylims=(0.995, 1.00))
scatter!(albedos[ind_lowAOD], ch4_mw2[ind_lowAOD]./n2o_mw2[ind_lowAOD], label="XCH₄/XN₂O; AOD = 0.05")
scatter!(albedos[ind_highAOD], ch4_mw2[ind_highAOD]./n2o_mw2[ind_highAOD], label="XCH₄/XN₂O; AOD = 0.18")
xlabel!("Albedo")
ylabel!("Retrieved/True Ratio")

scatter(albedos[ind_lowlowAOD], co2_mw1[ind_lowlowAOD]./n2o_mw1[ind_lowlowAOD], label="XCO₂/XN₂O; AOD = 0.01")
scatter!(albedos[ind_lowAOD], co2_mw1[ind_lowAOD]./n2o_mw1[ind_lowAOD], label="XCO₂/XN₂O; AOD = 0.05")
scatter!(albedos[ind_highAOD], co2_mw1[ind_highAOD]./n2o_mw1[ind_highAOD], label="XCO₂/XN₂O; AOD = 0.18")
xlabel!("Albedo")
ylabel!("Retrieved/True Ratio")

scatter(albedos[ind_lowlowAOD], n2o_mw1[ind_lowlowAOD]./n2o_mw2[ind_lowlowAOD], label="XCO₂/XN₂O; AOD = 0.01")
scatter!(albedos[ind_lowAOD], n2o_mw1[ind_lowAOD]./n2o_mw2[ind_lowAOD], label="XCO₂/XN₂O; AOD = 0.05")
scatter!(albedos[ind_highAOD], n2o_mw1[ind_highAOD]./n2o_mw2[ind_highAOD], label="XCO₂/XN₂O; AOD = 0.18")
xlabel!("Albedo")
ylabel!("Retrieved/True Ratio")