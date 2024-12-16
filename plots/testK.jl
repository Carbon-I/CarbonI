using HDF5
using LinearAlgebra
using InstrumentOperator
using DelimitedFiles
using Interpolations
using Statistics
FT = Float32

# Run using: /net/fluo/data2/software/Julia/julia-1.9.3/bin/julia --project=/home/sanghavi/code/github/OCORaman/OCORaman/ /home/sanghavi/code/github/OCORaman/OCORaman/src/OCOPlots/gain.jl &
#RetrievalResults
#vector_altitude_levels
#  vector_altitude_levels_apriori
#  vector_pressure_levels
#  vector_pressure_levels_apriori
#  vector_pressure_levels_met
function print_h5_structure(group, level=0)
    indent = "  " ^ level  # Indentation for nested levels
    for key in keys(group)
        println("$indent$key")
        if typeof(group[key]) == HDF5.Group  # Check if the object is a group
            print_h5_structure(group[key], level + 1)  # Recursively explore the group
        end
    end
end

#Constants for later OCO A-band noise computation
MaxMS = [7.0e20, 2.45e20, 1.25e20] # photons/m^2/sr/μm/s
#c_bkg = [] #[0.0042]
#c_pht = [] #[0.0089]

L1bSc_files = [#"oco2_L1bScGL_37243a_210702_B11006r_220323023559.h5"  
    #"oco2_L1bScGL_42843a_220722_B11008r_220831173918.h5"  
    "oco2_L1bScND_37214a_210630_B11006r_220326081005.h5"  
    "oco2_L1bScND_42610a_220706_B11008r_220831171738.h5"
    #"oco2_L1bScGL_42453a_220625_B11008r_220810234014.h5"  
    #"oco2_L1bScGL_44539a_221115_B11008r_221213003743.h5"  
    "oco2_L1bScND_41856a_220515_B11008r_220804213809.h5"  
    "oco2_L1bScND_45442a_230116_B11008r_230209231005.h5"]

L2Agg_files = [#"Glint_37243_TallVeg.h5"  
        #"Glint_42843_Himalayas.h5"  
        "Nadir_37214_TallVeg.h5"  
        "Nadir_42610_Himalayas.h5"
        #"Glint_42453_Sahara.h5"  
        #"Glint_44539_Andes.h5"  
        "Nadir_41856_Sahara.h5"  
        "Nadir_45442_Andes.h5"]

Output_files = [#"out_Glint_TallVeg.h5"  
        #"out_Glint_Himalayas.h5"  
        "out_Nadir_TallVeg.h5"  
        "out_Nadir_Himalayas.h5"
        #"out_Glint_Sahara.h5"  
        #"out_Glint_Andes.h5"  
        "out_Nadir_Sahara.h5"  
        "out_Nadir_Andes.h5"]

fctr = 3
 
    path = "/home/sanghavi/data/OCO2_jacobians/"
    L1bSc_file = L1bSc_files[fctr]#"oco2_L1bScND_41856a_220515_B11008r_220804213809.h5"
    #SliceMeasurements/radiance_slice_o2
    h5open(path*L1bSc_file, "r") do f
        global MM  = read(f["FootprintGeometry/footprint_stokes_coefficients"])
        global l1_sid = read(f["SoundingGeometry/sounding_id"])
        global c_snr_coef = read(f["InstrumentHeader/snr_coef"])  
        #global R_O2A  = read(f["SoundingMeasurements/radiance_o2"])
        #global R_WCO2 = read(f["SoundingMeasurements/radiance_weak_co2"])
        #global R_SCO2 = read(f["SoundingMeasurements/radiance_strong_co2"])  
    end

    # Get all files in the directory
    all_files = readdir(path)
    ststr = replace(L1bSc_file[1:end-15], "L1bSc"=>"L2Met")
    # Filter the files based on the desired pattern
    matching_files = filter(file -> startswith(file, ststr), all_files)
    L2Met_file = path*matching_files[1]
    @show L1bSc_file
    @show L2Met_file

    # L2 aggregate input
    #h5open("/home/sanghavi/data/OCO2_jacobians/Nadir_41856_Sahara.h5", "r") do f
    h5open("/home/sanghavi/data/OCO2_jacobians/"*L2Agg_files[fctr], "r") do f
        #println("HDF5 File Structure:")
        #print_h5_structure(f)

        global l2_sid = read(f["RetrievalHeader/sounding_id"])

        global K   = read(f["RetrievalResults/jacobian"])
        global Sₚ  = read(f["RetrievalResults/aposteriori_covariance_matrix"])
        global Sₐ  = read(f["RetrievalResults/apriori_covariance_matrix"])
        global A   = read(f["RetrievalResults/averaging_kernel_matrix"])

        global sza = read(f["RetrievalGeometry/retrieval_solar_zenith"])
        global vza = read(f["RetrievalGeometry/retrieval_zenith"])
        global saz = read(f["RetrievalGeometry/retrieval_solar_azimuth"])
        global vaz = read(f["RetrievalGeometry/retrieval_azimuth"])
        global lat = read(f["RetrievalGeometry/retrieval_latitude"])
        global lon = read(f["RetrievalGeometry/retrieval_longitude"])
        global PA  = read(f["RetrievalGeometry/retrieval_polarization_angle"])
        global alt = read(f["RetrievalGeometry/retrieval_altitude"])
        #global alt2  = read(f["RetrievalResults/vector_altitude_levels"])
        #global alt3  = read(f["RetrievalResults/vector_altitude_levels_apriori"])
        global xₐ  = read(f["RetrievedStateVector/state_vector_apriori"])
        global x   = read(f["RetrievedStateVector/state_vector_result"])
        global x_names = read(f["RetrievedStateVector/state_vector_names"])
        global y   = read(f["SpectralParameters/measured_radiance"]) 
        global Sϵ  = read(f["SpectralParameters/measured_radiance_uncert"])
        global Fx  = read(f["SpectralParameters/modeled_radiance"])
        global λ   = read(f["SpectralParameters/wavelength"])
        global dry_sub_columns = read(f["RetrievalResults/retrieved_dry_air_column_layer_thickness"])

        global sample_idx = read(f["SpectralParameters/sample_indexes"])
        global n = size(K)[1] # Number of state vector elements
        global m = size(K)[2] # Length of measurement vector
        global N = size(K)[3] # Number of ground pixels
    end

    dictFile = "/home/cfranken/code/gitHub/InstrumentOperator.jl/json/oco2.yaml"
    oco = InstrumentOperator.load_L1(dictFile,path*L1bSc_file,L2Met_file);
    # Pick some bands as tuple (or just one)
    bands=(1,); #bands = (1,2,3);


    
    fidx=[collect(22:41), collect(48:49), collect(54:56), collect(63:64)]
    fidx=vcat(fidx...) # array of state vector indices relevant to the O2 A-band 

    #=for ipix=1:1
        Ki = K[:,:,ipix]' # mxn
        Ai = A[:,:, ipix] # nxn
        Gi = Ai * pinv(Ki) # nxm
    end=#

    # Get Raman spectra Δy for SZA, psurf, and albedo
    a=[]
    for i=0:20
        push!(a, acosd(i/20))  
    end
    r_sza=reverse(Int.(ceil.(a[8:21])))
    r_ρ = zeros(FT,21)
    r_ρ_str = []
    for iρ = 1:21
        r_ρ[iρ] = (iρ-1)*0.05
        push!(r_ρ_str, replace(string(round(r_ρ[iρ], digits=2)),"."=>"p"))
    end
    r_psurf=[1000, 750, 500]
    r_psurf=reverse(r_psurf)
    wo_sza = findall(x -> x<=70.0, sza)
    N_ = length(wo_sza)
    N  = length(sza)
    vρ₁ = zeros(N) #zeros(size(x)[2])
    vρ₂ = zeros(N) #zeros(size(x)[2])
    vρ₃ = zeros(N) #zeros(size(x)[2])
    dρ₁ = zeros(N) #zeros(size(x)[2])
    dρ₂ = zeros(N) #zeros(size(x)[2])
    dρ₃ = zeros(N) #zeros(size(x)[2])
    dx1 = zeros(n,N) #zeros(size(x)[1], size(x)[2])
    dx2 = zeros(size(fidx)[1],N) #zeros(size(fidx)[1], size(x)[2])
    dx3 = zeros(n-2,N) #zeros(size(x)[1]-2, size(x)[2])
    vλ₁ = zeros(1016,N) #zeros(1016, size(x)[2])
    vλ₂ = zeros(1016,N) #zeros(1016, size(x)[2])
    vλ₃ = zeros(1016,N) #zeros(1016, size(x)[2])
    vdy = zeros(1016,N) #zeros(1016, size(x)[2])
    dy_wSIF = zeros(3*1016,N) #zeros(1016, size(x)[2])
    dy_woSIF = zeros(3*1016,N) #zeros(1016, size(x)[2])
    vy1 = zeros(1016,N) 
    vy2 = zeros(1016,N) 
    vy3 = zeros(1016,N) 
    #Δp = zeros(10)
    #dx1 = zeros(size(x)[1], 10)
    #dx2 = zeros(size(fidx)[1], 10)
    # Compute Δx = G.Δy
    #wo_sahel = findall(x -> 20<=x<=25, lat)
    #for ipix=1:2000:size(x)[2]
    #for ictr=1:N_ #size(x)[2]
    ipix = 7854
    #for ictr=1:length(wo_sahel) #size(x)[2] #
        #ipix=ictr # wo_sahel[ictr]
        #ipix=wo_sza[ictr]
        #@show "starting", ipix, ictr
        # Getting wavelengths in each band
        wo1 = findall(x -> .7<x<.8, λ[:,ipix])
        wo2 = findall(x -> 1.5<x<1.8, λ[:,ipix])
        wo3 = findall(x -> 2<x, λ[:,ipix])

        λ₁ = (x[48,ipix] .+ sample_idx[wo1,ipix] .* x[49,ipix]) # λ[wo1]
        λ₂ = (x[50,ipix] .+ sample_idx[wo2,ipix] .* x[51,ipix]) # λ[wo2]
        λ₃ = (x[52,ipix] .+ sample_idx[wo3,ipix] .* x[53,ipix]) # λ[wo3]

        # Getting BRDFs in each band
        ν₀₁ = 1e7/770.
        ν₀₂ = 1e7/1615.
        ν₀₃ = 1e7/2060.

        ν₁  = 1e4./λ₁ .- ν₀₁
        ν₂  = 1e4./λ₂ .- ν₀₂  
        ν₃  = 1e4./λ₃ .- ν₀₃

        ρ₀  = 0.05
        Θ   = -0.1
        κ   = 0.75

        g = acos(-cosd(sza[ipix])*cosd(vza[ipix]) + sind(sza[ipix])*sind(vza[ipix])*cosd(vaz[ipix]-saz[ipix]))
        F = (1-Θ^2)/(1 + Θ^2 - 2Θ*cos(π-g))^1.5
        G = ((tand(sza[ipix]))^2 + (tand(vza[ipix]))^2 + tand(sza[ipix])*tand(vza[ipix])*cosd(vaz[ipix]-saz[ipix]))^0.5
        R = (1 - ρ₀)/(1 + G)
        F *= ρ₀*(cosd(sza[ipix])*cosd(vza[ipix])*(cosd(sza[ipix])+cosd(vza[ipix])))^(κ-1)   
        F *= (1 + R)

        ρ₁ = (x[39,ipix] .+ ν₁ * x[40,ipix] .+ ν₁.^2 * x[41,ipix])*F 
        ρ₂ = (x[42,ipix] .+ ν₂ * x[43,ipix] .+ ν₂.^2 * x[44,ipix])*F 
        ρ₃ = (x[35,ipix] .+ ν₃ * x[46,ipix] .+ ν₃.^2 * x[47,ipix])*F   
        
        vρ₁[ipix] = mean(ρ₁)
        vρ₂[ipix] = mean(ρ₂)
        vρ₃[ipix] = mean(ρ₃)

        #getting Stokes coefficients from l1b file
        id = l2_sid[ipix]
        wo = findall(x-> x==id, l1_sid)

        #wo1 = findall(x -> .7<x<.8, λ[:,ipix])
        #Sounding index
        GeoInd = collect(Tuple(wo[1]))
        # Get instrument slit data for that sounding:
        oco_sounding = InstrumentOperator.getMeasurement(oco, bands, (sample_idx[wo1,ipix],), GeoInd);

        MM₁ = MM[:,1,wo[1][1],wo[1][2]]
        c_pht = c_snr_coef[1,:,wo[1][1],:]
        c_bkg = c_snr_coef[2,:,wo[1][1],:]

        x_psurf = x[22,ipix]/100. # Pa -> hPa
        x_ρ     = π*ρ₁ # to convert from surface reflectance to something resembling a Lambertian albedo
        x_ρ[x_ρ.>1].=1 # this is a temporary patch
        x_sza   = sza[ipix]

        #=i_sza0 = findall(x-> x<=x_sza,r_sza)[end]
        str_sza0 = string(r_sza[i_sza0])
        i_sza1 = findall(x-> x>x_sza,r_sza)[1]
        str_sza1 = string(r_sza[i_sza1])=#
        i_sza0 = findall(x-> cosd(x)>=cosd(x_sza),r_sza)[end]
        str_sza0 = string(r_sza[i_sza0])
        i_sza1 = findall(x-> cosd(x)<=cosd(x_sza),r_sza)[1]
        str_sza1 = string(r_sza[i_sza1])

        i_psurf0 = findall(x-> x<=x_psurf,r_psurf)[end]
        str_psurf0 = string(r_psurf[i_psurf0])
        if(str_psurf0!="1000")
            i_psurf1 = findall(x-> x>x_psurf,r_psurf)[1]
            str_psurf1 = string(r_psurf[i_psurf1])
        else
            i_psurf1 = i_psurf0
            str_psurf1 = str_psurf0
        end

        i_ρ0 = findall(x-> x<=minimum(x_ρ),r_ρ)[end]
        str_ρ0 = replace(string(round(r_ρ[i_ρ0], digits=2)),"."=>"p")
        if(str_ρ0!="1p0")
            i_ρ1 = findall(x-> x>=maximum(x_ρ),r_ρ)[1]
            str_ρ1 = replace(string(round(r_ρ[i_ρ1], digits=2)),"."=>"p")
        else
            i_ρ1 = i_ρ0
            str_ρ1 = "1p0"
        end
        @show ipix, ':', x_psurf, minimum(x_ρ),'-',maximum(x_ρ), x_sza

        isurf = 3
        iρ = 1
        iA = 1
        sza_str  = string(r_sza[iA])
        alb_str  = r_ρ_str[iρ]
        psurf_str= string(r_psurf[isurf])
        fname0   = "/home/sanghavi/RamanSIFgrid/raylSIF_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_nors_ABO2.dat"
        specNoRS = readdlm(fname0) 
        specNoRS = specNoRS[end:-1:1,:]
        wl_ = 1e7./specNoRS[:,1]
        ind_to_elim = findall(x->x>1e-6, wl_[1:end-2]-wl_[3:end]) #indices of duplicate entries (to be eliminated)
        r_wl = wl_[ filter(x->!(x in ind_to_elim), eachindex(wl_)) ]
        Δλ₁ = 5.0 #nm
        ind_hires = findall(x->x>λ₁[1]*1000-Δλ₁ && x<λ₁[end]*1000+Δλ₁, r_wl); 
        rΔy₁ = zeros(length(λ₁), length(r_psurf[i_psurf0:i_psurf1]), length(r_ρ_str[i_ρ0:i_ρ1]), length(r_sza[i_sza0:i_sza1]))  
        oco_Δy₁ = zeros(length(λ₁))
        for isurf = i_psurf0:i_psurf1 # 1:1 #2:2 # 
            for iρ = i_ρ0:i_ρ1 #3 #1:15
                for iA = i_sza0:i_sza1
                    sza_str  = string(r_sza[iA])
                    alb_str  = r_ρ_str[iρ]
                    psurf_str= string(r_psurf[isurf])
                    iBand = 1 # O2 A-band

                    fname1_0   = "/home/sanghavi/data/RamanSIFgrid/rayl_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_rrs_ABO2.dat"
                    specRRS_  = readdlm(fname1_0)
                    specRRS_ = specRRS_[end:-1:1,:]
                    specRRS0 = zeros(length(r_wl), size(specRRS_,2))
                    for i=1:size(specRRS_,2)
                        specRRS0[:,i] = specRRS_[ filter(x->!(x in ind_to_elim), eachindex(specRRS_[:,i])), i]
                    end
                    
                    fname1_00   = "/home/sanghavi/data/RamanSIFgrid/rayl_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_nors_ABO2.dat"
                    specnoRS_  = readdlm(fname1_00)
                    specnoRS_ = specnoRS_[end:-1:1,:]
                    specnoRS0 = zeros(length(r_wl), size(specnoRS_,2))
                    for i=1:size(specnoRS_,2)
                        specnoRS0[:,i] = specnoRS_[ filter(x->!(x in ind_to_elim), eachindex(specnoRS_[:,i])), i]
                    end
                    #ΔR = (specRRS0[ind_hires, 2]+specRRS0[ind_hires,5])*MM₁[1] +
                    #    (specRRS0[ind_hires, 3]+specRRS0[ind_hires,6])*MM₁[2] +
                    #    (specRRS0[ind_hires, 4]+specRRS0[ind_hires,7])*MM₁[3]
                    ΔR = (specRRS0[ind_hires, 2]+specRRS0[ind_hires,5]-specnoRS0[ind_hires,2])*MM₁[1] -
                        (specRRS0[ind_hires, 3]+specRRS0[ind_hires,6]-specnoRS0[ind_hires,3])*MM₁[2] +
                        (specRRS0[ind_hires, 4]+specRRS0[ind_hires,7]-specnoRS0[ind_hires,4])*MM₁[3]
                    # first convolve
                    interp_I = LinearInterpolation(r_wl[ind_hires]*1e-3, ΔR);
                    res = 0.001e-3;
                    off = 0.5e-3
                    wl = oco_sounding.ils[iBand].ν_out[1]-off:res:oco_sounding.ils[iBand].ν_out[end]+off;
                    #@show wl[1],wl[end], λ_grid[1],λ_grid[end]
                    I_wl = interp_I(wl);
                    tmp1   = InstrumentOperator.conv_spectra(oco_sounding.ils[iBand], wl, I_wl)
                    interp_ΔR = LinearInterpolation(λ₁, tmp1)
                    # then resample
                    tmp2 = interp_ΔR(λ₁)
                    rΔy₁[:, isurf-i_psurf0+1, iρ-i_ρ0+1, iA-i_sza0+1] = tmp2 
                end
            end
        end 
        @show x_psurf, minimum(x_ρ),'-',maximum(x_ρ), x_sza

        ϵ = 1e-3
        p_grid = round(Float64(r_psurf[i_psurf0]), digits=1):250:round(Float64(r_psurf[i_psurf1]), digits=1)+ϵ 
        ρ_grid = round(Float64(r_ρ[i_ρ0]), digits=2):0.05:round(Float64(r_ρ[i_ρ1]), digits=2)+ϵ
        #r_ρ[i_ρ0]:0.05:r_ρ[i_ρ1]+ϵ
        
        #=μ_grid=zeros(1+i_sza1-i_sza0)
        ii=0
        for μctr = i_sza1:-1:i_sza0
            ii+=1
            μ_grid[ii] = cosd(r_sza[μctr])
            @show μ_grid[ii], r_sza[μctr], μctr, ii
        end=#
        μ_grid = cosd(r_sza[i_sza1]):0.05:cosd(r_sza[i_sza0])+0.025
        if !(μ_grid[1]<=cosd(x_sza)<=μ_grid[2])
            μ_grid = reverse(cosd(r_sza[i_sza0]):-0.05:cosd(r_sza[i_sza1])-0.025)
        end
        
        if i_psurf0==i_psurf1
            global rΔy₁ = repeat(rΔy₁, 1, 2, 1, 1)
            global p_grid = r_psurf[i_psurf0]:10:r_psurf[i_psurf0]+10
        end
        if i_ρ0==i_ρ1
            global rΔy₁ = repeat(rΔy₁, 1, 1, 2, 1)
            global ρ_grid = r_ρ[i_ρ0]:ϵ:r_ρ[i_ρ0]+1.5ϵ
        end
        if i_sza0==i_sza1
            global rΔy₁ = repeat(rΔy₁, 1, 1, 1, 2)
            global μ_grid = cosd(r_sza[i_sza0]):ϵ:cosd(r_sza[i_sza0])+1.5ϵ
        end
        #wl_Δy_sitp = Vector{Any}(lBand1, undef)
        @show p_grid
        @show ρ_grid
        @show μ_grid

        for lctr=1:length(λ₁)
            #@show 0, lctr
            w1_Δy_itp = interpolate(rΔy₁[lctr,:,:,end:-1:1], BSpline(Cubic(Interpolations.Line(OnGrid()))))    
            #@show 1, lctr
            w1_Δy_sitp = Interpolations.scale(w1_Δy_itp, p_grid, ρ_grid, μ_grid)
            #@show 2, lctr
            oco_Δy₁[lctr] = w1_Δy_sitp(x_psurf<1000 ? x_psurf : 1000, x_ρ[lctr]<=1 ? x_ρ[lctr] : 1, cosd(x_sza))         
            #@show 3, lctr
        end         
        #@show ipix 
        # conversion from W/sr/μm/m² to Photons/s/sr/μm/m² ("Ph sec^{-1} m^{-2} sr^{-1} um^{-1}")
        h_Pl  = 6.62607015e-34 # J.s
        c_l   = 299792458 # m/s
        #computing energy per photons
        E = (h_Pl*c_l)./(λ₁*1e-6) # J = W.s    

        wo0 = findall(x->x!=-999999.0,Sϵ[:,ipix])
        # full state and measurement vectors
        Ki = K[:,wo0,ipix]' # mxn
        Ai = A[:,:, ipix] # nxn
        Gi = Ai * pinv(Ki) # nxm
        gi = inv(Ki'*Diagonal(1 ./ Sϵ[wo0,ipix].^2)*Ki + inv(Sₐ[:,:,ipix]))*Ki'*Diagonal(1 ./ Sϵ[wo0,ipix].^2)
        mm  = size(Ki)[1]
        m₁ = length(oco_Δy₁)
        dx1[:,ipix] = gi * vcat(oco_Δy₁./E, zeros(mm-m₁))
        dy_wSIF[1:size(Ki,1),ipix] = Ki*dx1[:,ipix]
        #Δp[ipix] = dx1[22,ipix]

        #Compute measurement noise
        #=ocoSig = (y[wo0,ipix]) #in Photons/m^2/sr/um
        ocoSNR = vcat(
            sqrt.(100*ocoSig[wo1].^2 ./ (MaxMS[1]*(c_bkg[sample_idx[wo1],1].^2 * MaxMS[1]/100 .+ c_pht[sample_idx[wo1],1].^2 .* ocoSig[wo1]))),
            sqrt.(100*ocoSig[wo2].^2 ./ (MaxMS[2]*(c_bkg[sample_idx[wo2],2].^2 * MaxMS[2]/100 .+ c_pht[sample_idx[wo2],2].^2 .* ocoSig[wo2]))),
            sqrt.(100*ocoSig[wo3].^2 ./ (MaxMS[3]*(c_bkg[sample_idx[wo3],3].^2 * MaxMS[3]/100 .+ c_pht[sample_idx[wo3],3].^2 .* ocoSig[wo3]))))
        #oco_noise = zeros(length(ocoSig))
        oco_noise = (ocoSig./ocoSNR)=#
        

        # Only A-band related state parameters and measurements 
        Ki = K[fidx,wo1,ipix]' # mxn
        Ai = A[fidx,fidx, ipix] # nxn
        Gi = Ai * pinv(Ki) # nxm
        gi = inv(Ki'*Diagonal(1 ./ Sϵ[wo1,ipix].^2)*Ki + inv(Sₐ[fidx,fidx,ipix]))*Ki'*Diagonal(1 ./ Sϵ[wo1,ipix].^2)
        #m  = size(Ki)[1]
        #m₁ = length(oco_Δy₁)
        dx2[:,ipix] = gi * (oco_Δy₁./E)

        # no-SIF state and vector
        Ki = K[1:62,wo0,ipix]' # mxn
        Ai = A[1:62,1:62, ipix] # nxn
        Gi = Ai * pinv(Ki) # nxm
        gi = inv(Ki'*Diagonal(1 ./ Sϵ[wo0,ipix].^2)*Ki + inv(Sₐ[1:62,1:62,ipix]))*Ki'*Diagonal(1 ./ Sϵ[wo0,ipix].^2)
        mm  = size(Ki)[1]
        m₁ = length(oco_Δy₁)
        dx3[:,ipix] = gi * vcat(oco_Δy₁./E, zeros(mm-m₁))
        dy_woSIF[1:size(Ki,1),ipix] = Ki*dx3[:,ipix]
        #Δp[ipix] = dx3[22,ipix]

        dρ₁[ipix] = mean(dx1[39,ipix] .+ ν₁ * dx1[40,ipix] .+ ν₁.^2 * dx1[41,ipix])*F 
        dρ₂[ipix] = mean(dx1[42,ipix] .+ ν₂ * dx1[43,ipix] .+ ν₂.^2 * dx1[44,ipix])*F 
        dρ₃[ipix] = mean(dx1[35,ipix] .+ ν₃ * dx1[46,ipix] .+ ν₃.^2 * dx1[47,ipix])*F   
        
        vλ₁[1:length(λ₁),ipix] .= λ₁
        vλ₂[1:length(λ₂),ipix] .= λ₂
        vλ₃[1:length(λ₃),ipix] .= λ₃
        vdy[1:length(λ₁),ipix] .= oco_Δy₁./E
        vy1[1:length(λ₁),ipix] .= y[wo1,ipix]
        vy2[1:length(λ₂),ipix] .= y[wo2,ipix]
        vy3[1:length(λ₃),ipix] .= y[wo3,ipix]
        @show "ending", ipix, ictr
    end
    

