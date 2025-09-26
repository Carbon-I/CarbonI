using Plots
using DelimitedFiles, ImageFiltering
using LegendrePolynomials, Distributions, SpecialFunctions

const cSqrtLn2divSqrtPi  = 0.469718639319144059835
const cLn2               = 0.6931471805599
const cSqrtLn2           = 0.8325546111577
const cSqrt2Ln2          = 1.1774100225


ci = readdlm("data/citb_srf_shape.txt")[:,1]
xx = collect(-10:0.1:9.99) 
xP = range(-1,1, length(ci))
legDegree = 20
K = zeros(length(ci), legDegree+1)
for i in 0:legDegree
    P = Pl.(xP, i)
    K[:,i+1] = P
end
fitX = K \ log10.(ci)

Plots.plot(xx, log10.(ci), label="CIDU Measurement", linewidth = 2)
Plots.plot!(xx, K*fitX)
# Define an instrument (Convolution of 2 Box kernels with a Gaussian):
Δwl = 0.005
FWHM  = 1.15  # Full Width at Half Maximum in nm 
SSI  = 0.7   # Spectral Sampling Interval in nm
# Define the kernels:
# Spectral response function of the slit (2*SSI)
kern1 = CarbonI.box_kernel(1.8*SSI, Δwl)
# Spectral response function of the optics (Gaussian)
kern2 = CarbonI.gaussian_kernel(FWHM, Δwl)
# Spectral response function of the pixel (1*SSI)
kern3 = CarbonI.box_kernel(SSI, Δwl)
# Combine the kernels:
kernf = imfilter(imfilter(kern1, kern2), kern3)

x1 = (kern1.offsets[1]+1:abs(kern1.offsets[1])-1)*Δwl/SSI
x2 = (kern2.offsets[1]+1:abs(kern2.offsets[1])-1)*Δwl/SSI
x3 = (kern3.offsets[1]+1:abs(kern3.offsets[1])-1)*Δwl/SSI
xf = (kernf.offsets[1]+1:abs(kernf.offsets[1])-1)*Δwl/SSI


Plots.plot(xx, ci, label="CIDU Measurement", linewidth = 2, yscale=:log10)
#Plots.plot!(x1, kern1.parent/maximum(kern1.parent), label="Box 1.8SSI")
Plots.plot!(x2, kern2.parent/maximum(kern2.parent), label="Gaussian 0.35nm")
#Plots.plot!(x3, kern3.parent/maximum(kern3.parent), label="Box 0.7SSI")
#Plots.plot!(xf, kernf.parent/maximum(kernf.parent), label="Combined 2 boxes * Gaussian", linewidth = 2)

# Set up Distributions (in pixel units)
GaussSRF   = Normal(0.0, 0.65)
LorentzSRF = Cauchy(0.0, 0.025)
α = 0.25
a2 = imfilter((1-α) * pdf.(LorentzSRF, xx), α*pdf.(GaussSRF, xx))


Plots.plot(xx, ci, label="CIDU Measurement", linewidth = 2, yscale=:log10)
Plots.plot!(xx.+0.1, a2./maximum(a2))

Plots.plot!(xx, pdf.(GaussSRF, xx)/maximum(pdf.(GaussSRF, xx)), label="Gaussian", linewidth = 2)
Plots.plot!(xx, pdf.(LorentzSRF, xx)/maximum(pdf.(LorentzSRF, xx))*0.3, label="Lorentzian", linewidth = 2)
Plots.ylims!(1e-5, 1.1)

function voigt(x, γ_d, γ_l)
    y = sqrt(cLn2) * γ_l / γ_d
    #return cSqrtLn2divSqrtPi / γ_d * real(erfc(cSqrtLn2 / γ_d * x ) + im * y)
    z = (cSqrtLn2 / γ_d) * x + im * y
    return cSqrtLn2divSqrtPi / γ_d * real(erfcx(-im*z))
end


function sim_SRF(x; grid=xx)
    γ_d = FWHM/(2*cSqrt2Ln2)
    vo = voigt.(grid, x[1], x[2])
    return vo./maximum(vo)
end