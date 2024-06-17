using NCDatasets

function loadGeos(fp::String, fp2::string, fp3::string)
    ds = NCDataset(fp)
    a  = ds["hyam"][:];
    b  = ds["hybm"][:];
    ai = ds["hyai"][:];
    bi = ds["hybi"][:];
    lat  = ds["lat"][:];
    lon  = ds["lon"][:];
    time = ds["time"][:];
    ch4 = ds["SpeciesConc_CH4"][:];
    co = ds["SpeciesConc_CO"][:];
    co2 = ds["SpeciesConc_CO2"][:];
    h2o = ds["SpeciesConc_H2O"][:];
    c2h6 = ds["SpeciesConc_C2H6"][:];
    n2o = ds["SpeciesConc_N2O"][:];
    lev = a/1000 +b;
    levi = ai/1000 +bi;
    dp = abs.(diff(levi));
    dp /= sum(dp)
    dpTrop = deepcopy(dp)
    # Need to adjust this per tropopause height in GeosChem!
    dpTrop[35:end] .= 0;
    dpTrop /= sum(dpTrop);

    dsMet = NCDataset(fp2)
    pbl  = dsMet["Met_PBLH"][:];
    Tair = dsMet["Met_T"][:];
    tropoPause = dsMet["Met_TropHt"][:];
    tropoLev   = dsMet["Met_TropLev"][:];

    dsAero = NCDataset(fp3)
    aod_dust = dsAero["AODDust"][:];
    aod_BCPI = dsAero["AODHyg550nm_BCPI"][:];
    aod_OCPI = dsAero["AODHyg550nm_OCPI"][:];
    aod_SALA = dsAero["AODHyg550nm_SALA"][:];
    aod_SALC = dsAero["AODHyg550nm_SALC"][:];
    aod_SO4  = dsAero["AODHyg550nm_SO4"][:];

    psurf1 = dsMet["Met_PS1WET"][:];
    psurf2 = dsMet["Met_PS2WET"][:];
    psurf = (psurf1 + psurf2) / 2;



    close(dsMet); close(ds); close(dsAero)
end

function getColumnAverage(gas, dp)
    xgas = zeros(size(gas,1),size(gas,2));
    for i=1:size(gas,1)
        for j=1:size(gas,2)
            xgas[i,j] = gas[i,j,:,1]' * dp
        end
    end
    return xgas
end

function getTroposphericColumnAverage(gas, dp, tropoLev)
    xgas = zeros(size(gas,1),size(gas,2));
    for i=1:size(gas,1)
        for j=1:size(gas,2)
            ind = convert(Int, (ceil(tropoLev[i,j,1])))
            xgas[i,j] = gas[i,j,1:ind,1]' * dp[1:ind] / sum(dp[1:ind]);
        end
    end
    return xgas
end

xch4 = zeros(size(ch4,1),size(ch4,2));
for i=1:size(ch4,1)
    for j=1:size(ch4,2)
        xch4[i,j] = ch4[i,j,:,1]' * dp
    end
end

heatmap(lon, lat, ch4[:,:,1,1]'*1e9)

for 
