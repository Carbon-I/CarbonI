using SpecialPolynomials

function forward_model_x(𝐱::AbstractArray{FT} ;sun = solarIrr, instrument=lociBox, sza=sza, vza = 0.0, profile=profile,σ_matrix=σ_matrix, wl=wl) where {FT}
    dims = size(σ_matrix)
	@show dims
    vmrs = reshape(𝐱[1:(dims[2]*dims[3])],(dims[2],dims[3]) )
    poly = Legendre(𝐱[dims[2]*dims[3]+1:end])
    #@show size(vmrs)
    # Air Mass Factor
    AMF = 1/cosd(sza);# + 1/cosd(vza)
    
    # Total sum of τ
    ∑τ = zeros(FT,size(σ_matrix,1))
	#@show size(vmrs,2)
    for i=1:size(vmrs,2)
     #   @show i, vmrs[end,i]
         ∑τ[:] += sum(σ_matrix[:,:,i] .* (vmrs[:,i] .* profile.vcd_dry)', dims=2)
         #@show i, maximum(∑τ)
    end
    # Transmission without Tsolar
    T = sun .* reverse(exp.(-AMF * ∑τ))
	#@show ∑τ
    T_conv = CarbonI.conv_spectra(instrument, wl, T)
    L = T_conv;
    # x-axis for polynomial [-1,1], enables legendre later:
    x_poly = CarbonI.rescale_x(instrument.ν_out)
   return L .* poly.(x_poly) 
end

function forward_model_x_(𝐱::AbstractArray{FT} ;sun = solarIrr,reflectance=refl, instrument=lociBox, sza=sza, vza=0.0, profile=profile,σ_matrix=σ_matrix, wl=wl) where {FT}
    dims = size(σ_matrix)
	# @show dims
    vmrs = reshape(𝐱[1:(dims[2]*dims[3])],(dims[2],dims[3]) )
    poly = Polynomial(𝐱[dims[2]*dims[3]+1:end])
    #@show size(vmrs)
    # Air Mass Factor
    AMF = 1/cosd(sza) + 1/cosd(vza)
    
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
   return L
end