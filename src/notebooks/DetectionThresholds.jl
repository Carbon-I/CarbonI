### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ e235e45c-75a2-446b-90a6-37813a972264
using Unitful


# ╔═╡ ed68a462-7c27-4885-b786-1fc670174d7f
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 2500px;
    	padding-left: max(160px, 20%);
    	padding-right: max(160px, 20%);
	}
</style>
"""

# ╔═╡ 71e36ef5-d026-4b43-8f1d-f2deb21c5342
md"""
# Determining precision requirements for flux thresholds

### In support of the Carbon-I STM and Mission Requirements

---
"""

# ╔═╡ 67bc318a-fec6-11ef-0282-7f4a8b27f56c
md"""
Following Jacob et al.(2016), the minimum detection limit $Q_{\min}$ is expressed as:

$$Q_{\min} = 
\underbrace{
\frac{M_{\mathrm{gas}}}{M_{\mathrm{DryAir}}}
\,\frac{\sigma_{\mathrm{gas}}\,p}{g}
}_{\sigma^*_{\mathrm{gas}}\ \mathrm{(amount/area)}}
\,\times\,
q\,U\,W, \tag{D-6}$$
where $U$ is wind speed, $W$ is pixel size, $q$ is a factor (2 for detection or 5 for quantification), and $\sigma^*_{\mathrm{gas}}$ (in $\mathrm{kg/m}^{2}$ or $\mathrm{mol/m}^2$) is the single-measurement precision for the column mass, i.e. the error in $\Omega$.

To get rid of mass conversions, we can express the same equation in units of moles as well, i.e. $\frac{M_{\mathrm{gas}}}{M_{\mathrm{DryAir}}}$ will be obsolete and we can express $\sigma^*_{\mathrm{gas}}$ in mol/m$^2$ by multiplying the total mol/m$^2$ of dry air with the precision in terms of mixing ratios. 
To keep it simple, we can use the total amount of molecules in the atmosphere (integrated from 0-to-of-atmosphere) as 41mol/m$^3$ (molar volume at standard condition) * 8000m (scale height) = 328kmol/m$^2$ of dry air.

"""

# ╔═╡ 7f766f90-9c7a-409f-af2d-b582cace7329
# How many moles in ine cubic meter at 25degrees and 1atm (using ideal gas law)
Nmol = 41u"mol/m^3"

# ╔═╡ 56a007ec-0854-4135-825b-3001a7098120
# Scale height of the atmosphere (typical)
H = 8000u"m"

# ╔═╡ ad9f5711-2010-403e-9853-80b7b6dad7ae
# Total moles per m²
molPerArea = Nmol * H

# ╔═╡ 56b6a75f-cc2f-4825-a41c-34805403b158
md"""
Here we basically invert the equation above and solve for $\sigma$ in ppb
"""

# ╔═╡ 3fd58c66-bd05-431a-aea9-2b94e2b416c4
# Here we basically invert the equation above and solve for $\sigma$ in ppb
function getRequiredPrecision(Qmin, molarMass, U, W, q; molPerArea=molPerArea)
	σ = Qmin / molPerArea / (U*W*q) / molarMass
	Unitful.uconvert(u"ppb", σ)
end
	

# ╔═╡ 8e46c1b6-36b2-4ea9-8d8e-07f41dfc635b
ReferenceWind = 2u"m/s"

# ╔═╡ b87e9502-e98b-45d4-87a6-8e56c9a927a2
GlobalGSD = 400u"m"

# ╔═╡ 4df80df5-e61a-4065-825c-a3b71c2dd0d3
TargetGSD = 50u"m"

# ╔═╡ efec3d94-0d85-4009-8d73-b854b045b489
# 2 for detection, 5 for quantification
q = 2

# ╔═╡ d0853eb8-f1a4-4efa-bc00-138e938366a3
begin
	# Requirements for CH4:
	mMassCH4        = 16u"g/mol"
	global_ch4_flux = 175u"kg/hr"
	target_ch4_flux = 65u"kg/hr"
end

# ╔═╡ 6877ee0e-edce-4382-9251-4ed902e2f1ca
begin
	mMassCO2 = 44u"g/mol"
	global_co2_flux = 100e3u"kg/hr"
	target_co2_flux = 50e3u"kg/hr"
end

# ╔═╡ 2a25e8b9-8c69-4d90-ad96-e1a53e2c8049
begin
	mMassCO  = 28u"g/mol"
	global_co_flux = 1750u"kg/hr"
	target_co_flux = 650u"kg/hr"
end

# ╔═╡ ab274634-fbf8-4748-a92e-5306da8edac5
md"""
---
### Precision Requirements for CH₄
"""

# ╔═╡ 2b2bcfaa-0aef-4d36-b917-da51a39b8bd8
sigma_global_ch4 = getRequiredPrecision(global_ch4_flux, mMassCH4, ReferenceWind, GlobalGSD, q)

# ╔═╡ 29c95218-b4f8-48f5-9642-18c50a7d0ace
sigma_target_ch4 = getRequiredPrecision(target_ch4_flux, mMassCH4, ReferenceWind, TargetGSD, q)

# ╔═╡ 59c13d06-5521-498d-8470-4c8872424566
md"""
---

### Precision Requirements for CO₂
"""

# ╔═╡ ce39c46d-cb13-49a5-adc2-d34c2cad8511
sigma_global_co2 = getRequiredPrecision(global_co2_flux, mMassCO2, ReferenceWind, GlobalGSD, q)

# ╔═╡ 7d56905f-f5b5-468a-970b-ccd89fddc715
sigma_target_co2 = getRequiredPrecision(target_co2_flux, mMassCO2, ReferenceWind, TargetGSD, q)

# ╔═╡ dd0fa2fc-a50d-4af4-aab2-be6754d08d91
md"""
---

### Precision Requirements for CO
"""

# ╔═╡ 948780c0-000f-4faf-b176-8809bfa0f352
sigma_global_co = getRequiredPrecision(global_co_flux, mMassCO, ReferenceWind, GlobalGSD, q)

# ╔═╡ 3f8efb40-aa68-490f-ac7c-9c2d0749b27a
sigma_target_co = getRequiredPrecision(target_co_flux, mMassCO, ReferenceWind, TargetGSD, q)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
Unitful = "~1.22.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.2"
manifest_format = "2.0"
project_hash = "d3faad2d9e5b4711ecf5d2410bef49c5521aff50"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "c0667a8e676c53d390a09dc6870b3d8d6650e2bf"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.22.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"
"""

# ╔═╡ Cell order:
# ╟─ed68a462-7c27-4885-b786-1fc670174d7f
# ╟─71e36ef5-d026-4b43-8f1d-f2deb21c5342
# ╟─67bc318a-fec6-11ef-0282-7f4a8b27f56c
# ╠═e235e45c-75a2-446b-90a6-37813a972264
# ╠═7f766f90-9c7a-409f-af2d-b582cace7329
# ╠═56a007ec-0854-4135-825b-3001a7098120
# ╠═ad9f5711-2010-403e-9853-80b7b6dad7ae
# ╟─56b6a75f-cc2f-4825-a41c-34805403b158
# ╠═3fd58c66-bd05-431a-aea9-2b94e2b416c4
# ╠═8e46c1b6-36b2-4ea9-8d8e-07f41dfc635b
# ╠═b87e9502-e98b-45d4-87a6-8e56c9a927a2
# ╠═4df80df5-e61a-4065-825c-a3b71c2dd0d3
# ╠═efec3d94-0d85-4009-8d73-b854b045b489
# ╠═d0853eb8-f1a4-4efa-bc00-138e938366a3
# ╠═6877ee0e-edce-4382-9251-4ed902e2f1ca
# ╠═2a25e8b9-8c69-4d90-ad96-e1a53e2c8049
# ╟─ab274634-fbf8-4748-a92e-5306da8edac5
# ╠═2b2bcfaa-0aef-4d36-b917-da51a39b8bd8
# ╠═29c95218-b4f8-48f5-9642-18c50a7d0ace
# ╟─59c13d06-5521-498d-8470-4c8872424566
# ╠═ce39c46d-cb13-49a5-adc2-d34c2cad8511
# ╠═7d56905f-f5b5-468a-970b-ccd89fddc715
# ╟─dd0fa2fc-a50d-4af4-aab2-be6754d08d91
# ╠═948780c0-000f-4faf-b176-8809bfa0f352
# ╠═3f8efb40-aa68-490f-ac7c-9c2d0749b27a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
