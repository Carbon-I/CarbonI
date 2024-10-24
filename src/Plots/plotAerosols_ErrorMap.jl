using Interpolations, JLD2,ColorSchemes, ImageFiltering, CairoMakie
@load "simulated_rads_all.jld2" R_conv_carbonI_dict
sorted_keys = sort(collect(keys(R_conv_carbonI_dict)));

paod    = [key[2] for key in sorted_keys];
aods    = [key[3] for key in sorted_keys];
albedos = [key[4] for key in sorted_keys];
szas    = [key[1] for key in sorted_keys];

# Find a low AOD:
index = 26
ind = findall(aods .== unique(aods)[index] .&& szas .== 20 .&& paod .== 850)
aod = unique(aods)[index]

rad = [maximum(R_conv_carbonI_dict[sorted_keys[i]]) for i in ind]

alb_grid = unique(albedos)
# Create Interpolator with albedo:
inter = LinearInterpolation(alb_grid, rad, extrapolation_bc=Interpolations.Flat())


using NCDatasets
f = Dataset("data/Sentinel2_LosAngeles_60m.nc")
la = f["data"][:]

@load "mw2_fits_all.jld2"

#alb_map = la[100:3000, 2490:-1:1000]
alb_map = la[200:1200, 1600:-1:600]

ind = findall(aods .== unique(aods)[index] .&& szas .== 20 .&& paod .== 850)


inter = LinearInterpolation(alb_grid, n2o_mw2[ind], extrapolation_bc=Interpolations.Flat())
inter2 = LinearInterpolation(alb_grid, ch4_mw2[ind], extrapolation_bc=Interpolations.Flat())
# Plot:
# Create the heatmap


cb = cgrad(:vik25,  15, categorical = true)
cbA = cgrad(:viridis,  15, categorical = true)
fig = Figure(size=(1800, 500), fontsize=20, font = :bold)
ax1 = Axis(fig[1, 1], title="Albedo map")

ci = CairoMakie.heatmap!(ax1,alb_map,colormap=cbA, colorrange=(0.0,0.4))
Colorbar(fig[1, 2], ci,label="Albedo")#, label="CH₄ error (ppb)")
ax2 = Axis(fig[1, 3], title="N₂O error (%)")
# Error needs to be a function of Albedo as well
er_n2o = CairoMakie.heatmap!(ax2, imfilter((inter.(alb_map) .+ 0.04*randn(size(alb_map)) .- 1.0) * 100,Kernel.box((11,11))) , colormap=cb, colorrange=(-1,1))
ax3 = Axis(fig[1, 4], title="CH₄ error (%)")
er_ch4 = CairoMakie.heatmap!(ax3, (inter2.(alb_map) .+ 0.003*randn(size(alb_map)) .-1.0) * 100, colormap=cb, colorrange=(-1,1))
Colorbar(fig[1, 5], er_n2o,label="Retrieval error (%)") #, label="CH₄ error (ppb)")
ax4 = Axis(fig[1, 6], title="CH₄/N₂O")
errors = (inter2.(alb_map)./inter.(alb_map) .-1)*100 .+ 0.22
errors[alb_map .<0.03] .= NaN
er_ratio = CairoMakie.heatmap!(ax4, errors, colormap=cb, colorrange=(-0.25,0.25))
Colorbar(fig[1, 7], er_ratio,label="Proxy Error (%)") #, label="CH₄ error (ppb)")

for ax in (ax1,ax2,ax3,ax4)
    hidespines!(ax)
    hidexdecorations!(ax)
    hideydecorations!(ax)
end
fig

@load "errors_albedo_all.jld2" resis

alb = 0.02:0.02:0.5
# Upscale to 60m already here (/sqrt(4))
n2o_precision = resis[1,:,5]/320/sqrt(4)
ch4_precision = resis[1,:,1]/1900/sqrt(4)

inter_n2o_prec = CubicSplineInterpolation(alb, n2o_precision, extrapolation_bc=Interpolations.Flat())
inter_ch4_prec = CubicSplineInterpolation(alb, ch4_precision, extrapolation_bc=Interpolations.Flat())

fig = Figure(size=(1600, 800), fontsize=20, font = :bold)
ax1 = Axis(fig[1, 1], title="Albedo map")

ci = CairoMakie.heatmap!(ax1,alb_map,colormap=cbA, colorrange=(0.0,0.4))
Colorbar(fig[1, 2], ci,label="Albedo")#, label="CH₄ error (ppb)")
ax2 = Axis(fig[1, 3], title="N₂O precision error (%)")
# Error needs to be a function of Albedo as well
er_n2o = CairoMakie.heatmap!(ax2, (inter_n2o_prec.(alb_map) ) * 100 , colormap=cbA, colorrange=(0,5)) # .* randn(size(alb_map))
ax3 = Axis(fig[1, 4], title="CH₄ precision error (%)")
er_ch4 = CairoMakie.heatmap!(ax3, (inter_ch4_prec.(alb_map) ) * 100, colormap=cbA, colorrange=(0,5)) # .* randn(size(alb_map))
Colorbar(fig[1, 5], er_n2o,label="Precision error (%)") #, label="CH₄ error (ppb)")
ax4 = Axis(fig[2, 1], title="CH₄/N₂O")
errors = (inter2.(alb_map)./inter.(alb_map) .-1)*100 .+ 0.22
errors[alb_map .<0.03] .= NaN
er_ratio = CairoMakie.heatmap!(ax4, errors, colormap=cb, colorrange=(-0.25,0.25))
Colorbar(fig[2, 2], er_ratio,label="Systematic Proxy Error (%)") #, label="CH₄ error (ppb)")

ax5 = Axis(fig[2, 3], title="N₂O aerosol systematic error (%)")
# Error needs to be a function of Albedo as well
er_n2oo = CairoMakie.heatmap!(ax5, (inter.(alb_map) .- 1.0) * 100 , colormap=cb, colorrange=(-2,2))
ax6 = Axis(fig[2, 4], title="CH₄ aerosol systematic error (%)")
er_ch4 = CairoMakie.heatmap!(ax6, (inter2.(alb_map).-1.0) * 100, colormap=cb, colorrange=(-2,2))
Colorbar(fig[2, 5], er_n2oo,label="Systematic error (%)") #, label="CH₄ error (ppb)")

for ax in (ax1,ax2,ax3,ax4, ax5, ax6)
    hidespines!(ax)
    hidexdecorations!(ax)
    hideydecorations!(ax)
end
fig





#Test Stuff
alb_data  = alb_map[100:500,100:500][:]
n2o_error = er[100:500,100:500][:]
bins = 0.05:0.01:0.25
m_error = zeros(length(bins)-1)
for i in eachindex(bins[1:end-1])
    ii = findall(bins[i].< alb_data[:] .< bins[i+1])
    if length(ii) > 10
        m_error[i] = mean(n2o_error[ii])
    else
        m_error[i] = NaN
    end
    
end