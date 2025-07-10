### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 8b1fbb9c-5cef-11f0-2b17-eda0e920ed3c
begin
        import Pkg
        Pkg.activate("../..");
end

# ╔═╡ 8db4f671-ec22-470b-bddf-2a846e4745c5
# We are using the CarbonI and vSmartMOM packages here
using CarbonI, vSmartMOM, Plots

# ╔═╡ 01c27b50-caf2-4d71-8a9e-c0e4d164cc35
using Interpolations, Statistics, NCDatasets, Unitful, InstrumentOperator, PlutoUI 

# ╔═╡ 46ef4760-83b5-4205-a525-693dfd09f4aa
html"""
<style>
        main {
                margin: 0 auto;
                max-width: 2000px;
        padding-left: max(160px, 10%);
        padding-right: max(160px, 10%);
        }
</style>
"""


# ╔═╡ 46e48450-8330-4ec6-bb26-5a6035a44090
plotlyjs()

# ╔═╡ aca33b65-c454-4d89-8dfd-1950245e58b6
begin
	sunFile = joinpath(dirname(pathof(CarbonI)),"../", "data/solar_irr.nc")
	DS = Dataset(sunFile)
	wlSol = 1e3*DS["wl"][:]
	solar_irr = 1e3*DS["solar_irr"][:] # convert to mW/m2/nm
	close(DS)
	sol  = CubicSplineInterpolation(range(wlSol[1],wlSol[end], length=length(wlSol)),solar_irr, extrapolation_bc=Interpolations.Flat());
end

# ╔═╡ 8c8535dd-9f83-4345-a8be-ea39c5f9bfad


# ╔═╡ f7deb652-1cf8-46c2-a289-4735955b3e88
wl = 2030:0.004:2390

# ╔═╡ ae9cfd12-a4af-4f7c-aaa5-341bb744a1a9
# Define an instrument:
cM, wl_ci = CarbonI.create_carbonI_conv_matrix_cbe(wl);

# ╔═╡ 6cf963be-b7c6-474c-b49b-3d97dea8d686
parameters = parameters_from_yaml("/home/cfranken/code/gitHub/CarbonI/src/yaml/carbon-i.yaml");

# ╔═╡ 05829a22-06d9-46a5-a72a-e59e25dd7f72
model = model_from_parameters(parameters);

# ╔═╡ cd2db74e-55fb-4cde-b045-0ebf9354da6f
# Run full multiple scattering radiative transfer code:
a = rt_run(model);

# ╔═╡ 071f7292-b3c3-4227-9ebc-4126a97fccb4
R = a[1][1,1,:];

# ╔═╡ 14f3bbd6-7cb2-4d9d-a62a-351657effc74
begin
	ν = parameters.spec_bands[1];
	Δν = mean(diff(ν));
	gridi = ν[1]:Δν:ν[end]+10eps();
end

# ╔═╡ 1ae1201f-860d-41d7-95d3-73b72e273618
rad_inter = CubicSplineInterpolation(gridi, R);

# ╔═╡ e162a431-aa60-476a-b97a-99f58ed3e042
F = cM*(rad_inter(1e7./wl).*sol(wl));

# ╔═╡ 97563e9b-e8ec-4727-b6b6-6d5c47ff5195
begin
	# Define the instrument:
	# ET  = 44.0u"ms"         # Exposure time
	SSI = (2*0.7)u"nm"      # Spectral resolution
	Pitch = 18.0u"μm"       # Pixel pitch
	FPA_QE = 0.85           # FPA quantum efficiency
	Bench_efficiency = 0.63 # Bench efficiency
	Fnumber = 2.2           # F-number
	#readout_noise = 100    # Readout noise
	#dark_current = 10e3u"1/s" # Dark current
end

# ╔═╡ 276bbd41-6b18-42d2-a630-85dd39b1fbbf
md"""
**Exposure time (ms):** $(@bind ET_ms Slider(1:100; default=44, show_value=true))
"""

# ╔═╡ 8f6674ce-c723-4be2-9130-9d3de4444ab1
ET = ET_ms  * u"ms";          # Exposure time with units

# ╔═╡ 50c8c69d-3480-45ad-9c5d-788815cec5f0
md"""
**Readout noise (e⁻):** $(@bind readout_noise Slider(50:150; default=100, show_value=true))
"""

# ╔═╡ 37b59a98-605d-4fab-b1d8-8f921b57e011
md"""
**Dark current (e⁻/s):** $(@bind dark_current_ Slider(1e3:100e3; default=5e3, show_value=true))
"""

# ╔═╡ 9ea7e2db-4c43-42f9-877b-9595abee54ab
dark_current = (dark_current_)u"1/s"

# ╔═╡ 2364456e-5674-469a-bc91-bf095b49abd1
begin
		ins  = InstrumentOperator.createGratingNoiseModel(ET, Pitch,FPA_QE, Bench_efficiency, Fnumber, SSI, (readout_noise), dark_current);    
		nesr = InstrumentOperator.noise_equivalent_radiance(ins, (wl_ci)u"nm", (F)u"mW/m^2/nm/sr");
		nesr_ = nesr./1u"mW/m^2/nm/sr"
		e = InstrumentOperator.photons_at_fpa(ins, (wl_ci)u"nm", (F)u"mW/m^2/nm/sr") .+ dark_current * ET;
end

# ╔═╡ de918053-64c4-47ef-9b64-ce7891da84ff
md"""
**Co-adding (1 for Target, 10 for Global):** $(@bind coadding Slider(1:1:20; default=10, show_value=true))
"""

# ╔═╡ bf350dde-8540-4172-8a0c-da771d4892c7
begin
	plot(wl_ci, F, label="Simulated Carbon-I measurement")
	xlabel!("Wavelength (nm)")
	ylabel!("Reflected Radiance (mW/m²/nm/sr)")
end

# ╔═╡ 2f870569-1fb2-4d59-aed6-cf8bd4dce85b
begin
    l  = @layout [a  b]

    p1 = plot(wl_ci, F ./ (nesr_/sqrt(coadding));
              label = "SNR",
              xlabel = "Wavelength (nm)",
              ylabel = "SNR")

    p2 = plot(wl_ci, e;
              label = "Electrons at FPA",
              xlabel = "Wavelength (nm)",
              ylabel = "Electrons at single FPA pixel")

    # wider canvas: 1200 px × 400 px
    plot(p1, p2; layout = l, size = (1200, 400))
end

# ╔═╡ Cell order:
# ╠═8b1fbb9c-5cef-11f0-2b17-eda0e920ed3c
# ╠═46ef4760-83b5-4205-a525-693dfd09f4aa
# ╠═8db4f671-ec22-470b-bddf-2a846e4745c5
# ╠═01c27b50-caf2-4d71-8a9e-c0e4d164cc35
# ╠═46e48450-8330-4ec6-bb26-5a6035a44090
# ╟─aca33b65-c454-4d89-8dfd-1950245e58b6
# ╠═8c8535dd-9f83-4345-a8be-ea39c5f9bfad
# ╠═f7deb652-1cf8-46c2-a289-4735955b3e88
# ╠═ae9cfd12-a4af-4f7c-aaa5-341bb744a1a9
# ╠═6cf963be-b7c6-474c-b49b-3d97dea8d686
# ╠═05829a22-06d9-46a5-a72a-e59e25dd7f72
# ╠═cd2db74e-55fb-4cde-b045-0ebf9354da6f
# ╠═071f7292-b3c3-4227-9ebc-4126a97fccb4
# ╟─14f3bbd6-7cb2-4d9d-a62a-351657effc74
# ╠═1ae1201f-860d-41d7-95d3-73b72e273618
# ╠═e162a431-aa60-476a-b97a-99f58ed3e042
# ╠═97563e9b-e8ec-4727-b6b6-6d5c47ff5195
# ╠═2364456e-5674-469a-bc91-bf095b49abd1
# ╟─9ea7e2db-4c43-42f9-877b-9595abee54ab
# ╟─8f6674ce-c723-4be2-9130-9d3de4444ab1
# ╟─276bbd41-6b18-42d2-a630-85dd39b1fbbf
# ╟─50c8c69d-3480-45ad-9c5d-788815cec5f0
# ╟─37b59a98-605d-4fab-b1d8-8f921b57e011
# ╟─de918053-64c4-47ef-9b64-ce7891da84ff
# ╟─bf350dde-8540-4172-8a0c-da771d4892c7
# ╟─2f870569-1fb2-4d59-aed6-cf8bd4dce85b
