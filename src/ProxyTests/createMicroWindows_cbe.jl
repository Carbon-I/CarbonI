
# Generic settings:
n_poly_2 = 5
n_poly_1 = 5
n_poly_3 = 8
n_poly_4 = 8

pbl_error = 10000.0
n_layers = 3
rel_errors_gases = 0.002
rel_n2o = 0.000001
n_iter = 5
# Something low here to account for our no noise scenario
Se_error = 0.00002


# Example Gas Array (we merge CO2_13 and CO2 now) ########################
gas_array1 = ["co2", "ch4", "n2o", "h2o","hdo"];
#gas_array = ["ch4", "co2", "n2o", "h2o", "co", "hdo", "c2h6"];
#gas_array = ["ch4", "co2", "n2o",  "co",  "c2h6"];
#gas_array = ["ch4", "n2o", "h2o", "co", "hdo", "c2h6"];
setupFile = "/home/cfranken/code/gitHub/CarbonI/src/yaml/carbon-i.yaml"
n_layers = 3
indLR1 = 8:150
#indLR = 287:410
#indLR = 7:410
# CLS = correlation length scale
cls        = Dict(gas => 550.0 for gas in gas_array1)
cls["h2o"] = 150.0
#cls["co2"] = 150.0
cls["hdo"] = 150.0
rel_errors = Dict(gas => rel_errors_gases for gas in gas_array1)
rel_errors["h2o"] = 0.5
rel_errors["hdo"] = 0.5
rel_errors["n2o"] = rel_n2o


############################################################################




a1 = createFitParams(indLR1, setupFile, n_layers, gas_array1)
xa1, Sa1, h_column1 = createBayesianConstraints(gas_array1, a1.gasProfiles,n_poly_1, a1.profile, 150.0, 800.0, cls, rel_errors, pbl_error)
ff = define_forward_model(a1)
errors = 1e10*ones(length(indLR1));
errors .= Se_error
Se1 = Diagonal(errors.^2);
w1 = findall(h_column1["n2o"].> 0.0)
Sa1[w1[end],w1[end]] = 100^2
#Sa[15,15] = 100^2
iP = length(xa1)-n_poly_1+1
ii_mw1 = define_inverse_model(ff,xa1,Sa1,Se1,length(indLR1),n_iter,iP);



# Example Gas Array (we merge CO2_13 and CO2 now) ########################
gas_array4 = ["co2", "ch4", "n2o", "h2o","hdo"];
#gas_array = ["ch4", "co2", "n2o", "h2o", "co", "hdo", "c2h6"];
#gas_array = ["ch4", "co2", "n2o",  "co",  "c2h6"];
#gas_array = ["ch4", "n2o", "h2o", "co", "hdo", "c2h6"];
setupFile = "/home/cfranken/code/gitHub/CarbonI/src/yaml/carbon-i.yaml"

indLR4 = 150:270
#indLR = 287:410
#indLR = 7:410
# CLS = correlation length scale
cls        = Dict(gas => 1550.0 for gas in gas_array4)
cls["h2o"] = 150.0
#cls["co2"] = 150.0
cls["hdo"] = 150.0
rel_errors = Dict(gas => rel_errors_gases for gas in gas_array4)
rel_errors["h2o"] = 0.5
rel_errors["hdo"] = 0.5
rel_errors["n2o"] =rel_n2o


############################################################################




a4 = createFitParams(indLR4, setupFile, n_layers, gas_array4)
xa4, Sa4, h_column4 = createBayesianConstraints(gas_array4, a4.gasProfiles,n_poly_4, a4.profile, 150.0, 800.0, cls, rel_errors, pbl_error)

ff = define_forward_model(a4)
errors = 1e10*ones(length(indLR4));
errors .= Se_error
Se4 = Diagonal(errors.^2);
w4 = findall(h_column4["n2o"].> 0.0)
Sa4[w4[end],w4[end]] = 100^2
iP = length(xa4)-n_poly_4+1
ii_mw4 = define_inverse_model(ff,xa4,Sa4,Se4,length(indLR4),n_iter,iP);


# Define window 2

# Example Gas Array (we merge CO2_13 and CO2 now) ########################
#gas_array = ["co2", "ch4", "n2o", "h2o","hdo"];
gas_array = ["ch4", "n2o", "h2o", "co","hdo", "c2h6"];
#gas_array = ["ch4", "co2", "n2o",  "co",  "c2h6"];
#gas_array = ["ch4", "n2o", "h2o", "co", "hdo", "c2h6"];
setupFile = "/home/cfranken/code/gitHub/CarbonI/src/yaml/carbon-i.yaml"

#indLR = 8:260
indLR2 = 270:399
#indLR2 = 300:400
#indLR = 7:410
cls        = Dict(gas => 1550.0 for gas in gas_array)
cls["h2o"] = 150.0
#cls["co2"] = 150.0
cls["hdo"] = 150.0
rel_errors = Dict(gas => rel_errors_gases for gas in gas_array)
rel_errors["h2o"] = 0.5
rel_errors["hdo"] = 0.5
rel_errors["n2o"] = rel_n2o


############################################################################



a = createFitParams(indLR2, setupFile, n_layers, gas_array)
xa2, Sa2, h_column2 = createBayesianConstraints(gas_array, a.gasProfiles,n_poly_2, a.profile, 150.0, 850.0, cls, rel_errors, pbl_error)
w = findall(h_column2["n2o"].> 0.0)
#Sa2[w[3],w[3]] = 10^2
ff = define_forward_model(a)
errors = 1e10*ones(length(indLR2));
errors .= Se_error
Se2 = Diagonal(errors.^2);
w2 = findall(h_column2["n2o"].> 0.0)
Sa2[w2[end],w2[end]] = 100^2
#Sa[15,15] = 100^2
iP = length(xa2)-n_poly_2+1
ii_mw2 = define_inverse_model(ff,xa2,Sa2,Se2,length(indLR2),n_iter,iP);



# Third window

# Example Gas Array (we merge CO2_13 and CO2 now) ########################
#gas_array = ["co2", "ch4", "n2o", "h2o","hdo"];
gas_array = ["ch4", "h2o", "co","hdo", "c2h6"];
#gas_array = ["ch4", "co2", "n2o",  "co",  "c2h6"];
#gas_array = ["ch4", "n2o", "h2o", "co", "hdo", "c2h6"];
setupFile = "/home/cfranken/code/gitHub/CarbonI/src/yaml/carbon-i.yaml"

#indLR = 8:260
indLR3 = 399:490
#indLR2 = 300:400
#indLR = 7:410
cls        = Dict(gas => 1550.0 for gas in gas_array)
cls["h2o"] = 150.0
#cls["co2"] = 150.0
cls["hdo"] = 150.0
rel_errors = Dict(gas => rel_errors_gases for gas in gas_array)
rel_errors["h2o"] = 0.5
rel_errors["hdo"] = 0.5
rel_errors["n2o"] = rel_n2o


############################################################################



a3 = createFitParams(indLR3, setupFile, n_layers, gas_array)
xa3, Sa3, h_column3 = createBayesianConstraints(gas_array, a3.gasProfiles,n_poly_3, a.profile, 150.0, 850.0, cls, rel_errors, pbl_error)
ff = define_forward_model(a3)
errors = 1e10*ones(length(indLR3));
errors .= Se_error
Se = Diagonal(errors.^2);
#Sa[15,15] = 100^2
iP = length(xa3)-n_poly_3+1
ii_mw3 = define_inverse_model(ff,xa3,Sa3,Se,length(indLR3),n_iter,iP);



