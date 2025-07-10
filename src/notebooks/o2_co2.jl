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

# ╔═╡ 0e4f9174-5495-488c-88a6-3495d0ebb007
begin
        import Pkg
        Pkg.activate("../..");
end

# ╔═╡ 4c0f0444-5d0d-11f0-36b8-a9cd579b6bf3
using Plots, DelimitedFiles, PlutoUI 

# ╔═╡ e8668fe0-0a0a-4381-9a7f-31b8658416af
using Polynomials          # polynomial fitting & evaluation

# ╔═╡ be169cd6-b788-44e4-857e-3dd9946484bb
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

# ╔═╡ 07401d16-a81e-4fae-b8c2-5c6cf9839285
gr()

# ╔═╡ 3b3cd743-b072-41ea-8b50-97223b46a15b
o2_cgo = readdlm("../../data/scripps/monthly_o2_cgo.csv", ',',skipstart=23);

# ╔═╡ fb899112-7c4e-4841-aa54-a72686145255
co2_cgo = readdlm("../../data/scripps/monthly_co2_cgo.csv", ',',skipstart=23);

# ╔═╡ b6fc2a1e-483b-4801-b5fc-fde9c68f8449
co2_mlo = readdlm("../../data/scripps/monthly_co2_mlo.csv", ',',skipstart=23);

# ╔═╡ 57e05f78-97b8-4335-b0d0-c9fcdf319d0a
o2_mlo = readdlm("../../data/scripps/monthly_o2_mlo.csv", ',',skipstart=23);

# ╔═╡ 943f25e8-92cd-4831-9db3-fb2ea2879337
# ╠═╡ disabled = true
#=╠═╡

  ╠═╡ =#

# ╔═╡ 1fa9da01-5bb0-4cb5-9c94-aa4174717839


# ╔═╡ 52f34139-49d3-45d0-b1b7-416e7523fbfb
# 1 ppm means 1 mole of CO₂ per 10⁶ moles of air.
# Using the molar mass ratio (CO₂ = 44 g mol⁻¹, mean air ≈ 28.97 g mol⁻¹) and the 5.15 × 10¹⁸ kg atmosphere gives 7.8 × 10¹⁵ kg = 7.8 Gt CO₂ per ppm.
ppm_per_gt = 1/7.8;

# ╔═╡ 886f006b-9d0f-46af-b0a8-6f874148cd39
function filter_data(A,B,C, D)
	# --- 2. Build a Boolean mask: true where neither value is –99.9 ------------
	mask = (A .!== -99.99) .& (B .!= -99.99) .& (C .!== -99.99) .& (D .!= -99.99)    # element-wise .!= and .&
	# --- 3. Get the matching indices -------------------------------------------
	idx = findall(mask)                      # vector of indices / CartesianIndex
end

# ╔═╡ efe41935-4cd8-485c-a6b8-5a7bbc36ff30
ind_all = filter_data(o2_cgo[:,5], co2_cgo[:,5], o2_mlo[:,5], co2_mlo[:,5]);

# ╔═╡ 8c0afa0e-e1f1-4ba1-9bba-27b5ecd94598
start = o2_cgo[ind_all[1],4]

# ╔═╡ 70279880-42e6-4227-94af-c09d442feede
n_years = o2_cgo[ind_all[end],4]-start

# ╔═╡ b6af9bcb-6423-4e96-9855-17986ca92136
subset = 1:350;

# ╔═╡ 2938bdb4-cc46-4760-869a-82c7e014d3e2
years = o2_cgo[ind_all[subset],4].-start;

# ╔═╡ 031ee78c-d887-48c5-81c2-58ae8df8d028
x = (co2_cgo[ind_all[subset],6] + co2_mlo[ind_all[subset],6])/2;

# ╔═╡ 2a8e487b-05c7-4294-995c-702d00703c6f
y = 0.5*(o2_cgo[ind_all[subset],6] + o2_mlo[ind_all[subset],6])  * 0.21;

# ╔═╡ 1f287372-b3b0-47e8-9dc4-2930c04043c1
p_lin  = fit(x, y, 1)      # degree-1 (linear)   =>  Poly([b, m])

# ╔═╡ 63196254-becb-4f69-82ab-c35f8c231c58
md"""
**Gigaton CO2 per year (on average):** $(@bind gt_per_year Slider(29:0.1:35; default=32, show_value=true))
"""

# ╔═╡ c2ae14d7-bebb-4288-9e89-d093e716798a
d_CO2_ppm = ppm_per_gt * gt_per_year *  years .+ x[1];

# ╔═╡ 9da99816-c5e8-4060-8d5a-71b36cfe0b28
md"""
**O₂ per CO₂ for combustion (mol/mol) :** $(@bind o2_per_co2 Slider(-1.45:0.01:-1.25; default=-1.33, show_value=true))
"""

# ╔═╡ 8ce0db19-4917-4fa7-a5a9-97444f618240
d_O2_ppm = ppm_per_gt * gt_per_year * years * o2_per_co2 .+ y[1];

# ╔═╡ bd9911ea-2df8-4d32-a3a5-124b65e6a6cc
md"""
**O₂ per CO₂ for photosynthesis (mol/mol) :** $(@bind o2_per_co2_photo Slider(-1.2:0.01:-1.0; default=-1.1, show_value=true))
"""

# ╔═╡ 5fa8b86d-b608-4935-a735-818c660b2e17
function separate_uptake(co2_end_measured, o2_end_measured, co2_end_fossil, o2_end_fossil, co2_start_fossil )
	dCO2 = co2_end_fossil - co2_end_measured
	dO2  = o2_end_fossil - o2_end_measured
	slope = dCO2/dO2
	land_uptake  = dO2 / o2_per_co2_photo
	ocean_uptake = dCO2 - land_uptake
	airborneFraction = (co2_end_fossil - co2_end_measured)/(co2_end_fossil-co2_start_fossil)
	return land_uptake, ocean_uptake, airborneFraction
end

# ╔═╡ 322f1162-59ae-4df5-b2cd-4877be7cde0f
land_sink, ocean_sink, af = separate_uptake(x[end], y[end], d_CO2_ppm[end], d_O2_ppm[end], d_CO2_ppm[1])

# ╔═╡ e929abe0-fbd2-4feb-ac0d-cb99910a639b


# ╔═╡ e3925d18-6f4f-4db0-bf2c-4304eba7fc74
begin
	plt = scatter(x, y;
	              label      = "Scripps data",
	              xlabel     = "CO₂ (ppm)",
	              ylabel     = "O₂ anomaly (per meg)",
	              legend     = :topright,
	              markerstrokecolor = :black,
	              markercolor       = :black)
	
	
	
	################################################################################
	# 2 ▪  “Fossil-fuel burning” arrow (black dashed guideline + red arrow)
	################################################################################
	lower_end = d_O2_ppm[end]-15
	
	# thick red arrow that stops at 2022 point
	plot!(d_CO2_ppm, d_O2_ppm;
	      lw    = 4,
	      color = :red,
	      label = "Fossil-fuel burning",
	      arrow = (:closed, 30))                 # arrowhead size ~10 px
	
	################################################################################
	# 3 ▪  Land-sink arrow  (green)
	################################################################################
	x_land_end = x[end] + land_sink
	y_land_end = y[end] + land_sink * o2_per_co2_photo   # photosyn. ratio
	
	plot!(reverse([x[end], x_land_end]), reverse([y[end], y_land_end]);
	      lw    = 4,
	      color = :forestgreen,
	      label = "Land uptake",
	      arrow = (:closed, 30))
	
	################################################################################
	# 4 ▪  Ocean-sink arrow  (blue, horizontal)
	################################################################################
	x_ocean_end = x_land_end + ocean_sink
	y_ocean_end = y_land_end         # flat (no O₂ change)
	
	plot!([x_ocean_end,x_land_end], [y_ocean_end,y_land_end];
	      lw    = 4,
	      color = :blue,
	      label = "Ocean uptake",
	      arrow = (:closed, 30));

	########################################################################
#  Values you already have
########################################################################
Δ_ocean = ocean_sink                         # the number you want to show
Δ_land  = land_sink
Δ_ff    = x[end] - x[1]                     # fossil-fuel ΔCO₂ you plotted
Δ_all   = 	Δ_ocean + Δ_land + Δ_ff
########################################################################
#  1 ▪  Ocean–sink arrow (blue, two directions)
########################################################################
plot!([x_ocean_end, x_land_end],
      [y_ocean_end-10, y_land_end-10];
      lw = 3, c = :blue, arrow = true, legend = false)
plot!(reverse([x_ocean_end, x_land_end]),
      [y_ocean_end-10, y_land_end-10];
      lw = 3, c = :blue, arrow = true, legend = false)

# label in the middle
x_mid_ocean = (x_ocean_end + x_land_end) / 2
y_mid_ocean = y_land_end - 10
annotate!(x_mid_ocean, y_mid_ocean + 3,             # small upward offset
          text("Ocean sink $(round(Δ_ocean/Δ_all*100; digits = 1))%", :center, 9, :blue))

########################################################################
#  2 ▪  Land-sink arrow (green)
########################################################################
plot!([x[end], x_land_end],
      [y_ocean_end-10, y_land_end-10];
      lw = 3, c = :green, arrow = true, legend = false)
plot!(reverse([x[end], x_land_end]),
      [y_ocean_end-10, y_land_end-10];
      lw = 3, c = :green, arrow = true, legend = false)

x_mid_land = (x[end] + x_land_end) / 2
y_mid_land = y_land_end - 10
annotate!(x_mid_land, y_mid_land + 3,
          text("Land Sink $(round(Δ_land/Δ_all*100; digits = 1))%", :center, 9, :forestgreen))

########################################################################
#  3 ▪  Fossil-fuel arrow (black)
########################################################################
plot!([x[1], x[end]],
      [y_ocean_end-10, y_land_end-10];
      lw = 3, c = :black, arrow = true, legend = false)
plot!(reverse([x[1], x[end]]),
      [y_ocean_end-10, y_land_end-10];
      lw = 3, c = :black, arrow = true, legend = false)

x_mid_ff = (x[1] + x[end]) / 2
y_mid_ff = y_land_end - 10
annotate!(x_mid_ff, y_mid_ff + 3,
          text("Airborne Fraction $(round(Δ_ff/Δ_all*100; digits = 1))%", :center, 9, :black))
	
	################################################################################
	# 5 ▪  Helper grid lines (optional)
	################################################################################
	# vertical dashed line through first CO₂ point
	vline!([x[1]]; lc = :gray, ls = :dash, lw=1.5, label = "")
	# horizontal dashed line through last O₂ point
	plot!([x[end]; x[end]], [lower_end; y[end]]; lc = :gray, ls = :dash, lw=1.5, label = "")
	plot!([d_CO2_ppm[end]; d_CO2_ppm[end]], [lower_end; d_O2_ppm[end]]; lc = :gray, ls = :dash, lw=1.5, label = "")
	plot!([x_land_end; x_land_end], [lower_end; d_O2_ppm[end]]; lc = :gray, ls = :dash, lw=1.5, label = "")
	
	################################################################################
	# 6 ▪  Final tweaks
	################################################################################
	#xlims!(minimum(x)-2,  x_ocean_end+5)
	#ylims!(minimum(y)-5,  maximum(y)+2)
	ylims!(lower_end, -10)
	plot!(size = (850, 750))          # square canvas like the original
end

# ╔═╡ Cell order:
# ╟─0e4f9174-5495-488c-88a6-3495d0ebb007
# ╟─be169cd6-b788-44e4-857e-3dd9946484bb
# ╠═4c0f0444-5d0d-11f0-36b8-a9cd579b6bf3
# ╠═e8668fe0-0a0a-4381-9a7f-31b8658416af
# ╠═07401d16-a81e-4fae-b8c2-5c6cf9839285
# ╠═3b3cd743-b072-41ea-8b50-97223b46a15b
# ╠═fb899112-7c4e-4841-aa54-a72686145255
# ╠═b6fc2a1e-483b-4801-b5fc-fde9c68f8449
# ╠═57e05f78-97b8-4335-b0d0-c9fcdf319d0a
# ╠═70279880-42e6-4227-94af-c09d442feede
# ╠═943f25e8-92cd-4831-9db3-fb2ea2879337
# ╠═8c0afa0e-e1f1-4ba1-9bba-27b5ecd94598
# ╠═1fa9da01-5bb0-4cb5-9c94-aa4174717839
# ╠═2938bdb4-cc46-4760-869a-82c7e014d3e2
# ╠═52f34139-49d3-45d0-b1b7-416e7523fbfb
# ╠═c2ae14d7-bebb-4288-9e89-d093e716798a
# ╠═8ce0db19-4917-4fa7-a5a9-97444f618240
# ╟─886f006b-9d0f-46af-b0a8-6f874148cd39
# ╠═efe41935-4cd8-485c-a6b8-5a7bbc36ff30
# ╠═b6af9bcb-6423-4e96-9855-17986ca92136
# ╠═031ee78c-d887-48c5-81c2-58ae8df8d028
# ╠═2a8e487b-05c7-4294-995c-702d00703c6f
# ╠═1f287372-b3b0-47e8-9dc4-2930c04043c1
# ╠═322f1162-59ae-4df5-b2cd-4877be7cde0f
# ╟─5fa8b86d-b608-4935-a735-818c660b2e17
# ╟─63196254-becb-4f69-82ab-c35f8c231c58
# ╟─9da99816-c5e8-4060-8d5a-71b36cfe0b28
# ╠═bd9911ea-2df8-4d32-a3a5-124b65e6a6cc
# ╟─e929abe0-fbd2-4feb-ac0d-cb99910a639b
# ╟─e3925d18-6f4f-4db0-bf2c-4304eba7fc74
