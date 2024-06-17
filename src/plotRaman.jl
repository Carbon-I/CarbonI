using Glob, DelimitedFiles, NCDatasets

files = glob("raylSIF*psurf1000hpa*nors_ABO2*.dat", "/home/sanghavi/RamanSIFgrid/");

n = length(files);
cRad = zeros(n);
inFill = zeros(n);
sza_values = zeros(n)
alb_values = zeros(n)


for (i, file) in enumerate(files)
    @show file
    specNoRS = readdlm(file);
    specRRS =  readdlm(replace(file, "nors" => "rrs"));
    cRad[i] = specRRS[5000,2];
    inFill[i] = specRRS[5000,5];
    (sza_values[i], alb_values[i]) = extract_values(file)
end

function extract_values(filename)
    # Regular expression to find the SZA value
    sza_pattern = r"sza(\d+)"
    sza_match = match(sza_pattern, filename)
    sza_value = sza_match !== nothing ? parse(Float64, sza_match.captures[1]) : NaN
    
    # Regular expression to find the albedo value (with 'p' as the decimal point)
    alb_pattern = r"alb(\d+)p(\d+)"
    alb_match = match(alb_pattern, filename)
    if alb_match !== nothing
        # Construct the float from parts before and after the 'p'
        alb_value = parse(Float64, alb_match.captures[1] * "." * alb_match.captures[2])
    else
        alb_value = NaN
    end
    
    return sza_value, alb_value
end

sza_values = zeros(n)
alb_values = zeros(n)

albs = reshape(alb_values,(15,14))[:,1]
szas = reshape(sza_values,(15,14))[1,:]

sorted_idx_x = sortperm(szas)
sorted_idx_y = sortperm(albs)

inFills = reshape(inFill,(15,14))
Rads = reshape(cRad,(15,14))

preFac = 1e7/755^2
contourf(szas[sorted_idx_x],albs[sorted_idx_y],preFac*inFills[sorted_idx_y,sorted_idx_x],levels=10)
xlabel!("SZA")
ylabel!("Albedo")
title!("Pseudo-SIF through Raman (mW/m²/sr/nn)")


nc = "/net/fluo/data3/DATASERVER/satellite/OCO2/L2_Lite_SIF/B11012Ar/original/2022/oco2_LtSIF_220525_B11012Ar_220822215313s.nc4"

ds = Dataset(nc)

rad_755 = ds.group["Science"]["continuum_radiance_757nm"]
relSIF_raw_755 = ds.group["Science"]["SIF_Unadjusted_Relative_757nm"]
relSIF_755 = ds.group["Science"]["SIF_Relative_757nm"]

bins = ds.group["Offset"]["signal_histogram_bins"][:]

corrFactor = zeros(length(bins))    
for i in eachindex(bins)
    ii = findall(x->x>bins[i]-3 && x<bins[i]+3, rad_755);
    if length(ii) > 0
        corrFactor[i] = mean(relSIF_raw_755[ii]-relSIF_755[ii]) 
    end
end

corrFactorRaman = zeros(length(bins))    
for i in eachindex(bins)
    ii = findall(x->x>bins[i]-3 && x<bins[i]+3, preFac*contRad[1,:,:]);
    if length(ii) > 0
        corrFactorRaman[i] = mean(addRaman[1,:,:][ii]./contRad[1,:,:][ii]) 
    end
end

corri = ds.group["Offset"]["SIF_Mean_757nm"][:]

L1File   = "/net/fluo/data3/data/FluoData1/group/oco2/L1bSc/oco2_L1bScND_26780a_190715_B10003r_200429212407.h5"
metFile  = "/net/fluo/data3/data/FluoData1/group/oco2/L2Met/oco2_L2MetND_26780a_190715_B10003r_200429212406.h5"
dictFile = "/home/cfranken/code/gitHub/InstrumentOperator.jl/json/oco2.yaml"

oco = InstrumentOperator.load_L1(dictFile,L1File, metFile);

# Pick some bands as tuple (or just one)
bands = (1,);
#bands = (1,3);
# Indices within that band:
indices = (2:1016,);
#indices = (92:885,50:916);
# Geo Index (footprint,sounding):
GeoInd = [5,5000];

# Get data for that sounding:
oco_sounding = InstrumentOperator.getMeasurement(oco, bands, indices, GeoInd);

