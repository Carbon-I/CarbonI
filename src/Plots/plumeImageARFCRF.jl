using NCDatasets, ImageFiltering
using StatsBase, LaTeXStrings, CairoMakie

les = Dataset("/net/fluo/data1/ftp/XYZT_LES_/Data_Train_Val_Test/original_train_set.h5")
#les = Dataset("/net/fluo/data2/groupMembers/siraput/Data_NewSatRetrieval/LES_P_T_W/C_conc_27777.h5")
#les = Dataset("/net/fluo/data2/groupMembers/siraput/Data_NewSatRetrieval/LES_P_T_W/run10_C/run10_C_26902.h5")
C3D = les.group["train_set"]["C_2D"][:,:,15,2];
C_sum = sum(C3D, dims=1)[1,:,:];

# From Siraput, unit conversion into ppm-m (in /net/fluo/data1/ftp/XYZT_LES_/code/helper_functions.py)
c2ppm = (500.0/3.57e-4)
c2mol = c2ppm*44/1e6
C_sum = C3D*10*c2mol # Use 200kg/hr (ref is about 5g/s


global_gsd = 330
target_gsd = 35.
target_alongTrack = 30.
fwhm_to_sigma = 2.35482004503
target_arf = 2.2/fwhm_to_sigma
target_crf = 1.3/fwhm_to_sigma
fac_global  = convert(Int,global_gsd/5)+1
kernel_global = Kernel.box((fac_global,fac_global))
kernel_target = Kernel.gaussian(((target_gsd/5)*target_crf,target_alongTrack/5*target_arf))

dim_x = (1:size(C_sum,1))*5 .-750
dim_y = (1:size(C_sum,1))*5 .-750

scene_global = imfilter(C_sum,kernel_global);
scene_target = imfilter(C_sum,kernel_target);


refColumn  = 4e19/6e23*100^2
# Using 4ppb here as baseline (15% albedo)
noise_target  = 4 * sqrt(10) / 1900 * refColumn
noise_global  = 4 * sqrt(10) / 1900 * refColumn


f = Figure(resolution=(450,400), title="", fontsize=16,backgroundcolor = :transparent)
ax1 = Axis(f[1,1],xlabel = "meters",ylabel = "meters", aspect = DataAspect(),title="Carbon-I Target Mode (200kg/hr)") 
sample=convert(Int,35/5)
hm = CairoMakie.heatmap!(ax1, dim_x[1:sample:end],dim_y[1:sample:end],  scene_target[1:sample:end,1:sample:end]./noise_target, colormap=:OrRd_9)
#hideydecorations!(ax1)
#hidexdecorations!(ax1)
#CairoMakie.ylims!(1e-5,1)
CairoMakie.xlims!(-200,700)
CairoMakie.ylims!(-200,700)
cb = Colorbar(f[1,2], hm,  label = L"\text{\Delta \Omega(CH_4)/1\sigma}",  labelpadding = 10, ticklabelpad = 2, width = 10)
f
save("plots/PixelEnhancements30M.pdf",f)


f = Figure(resolution=(450,400), title="Vertical Optical Depth of Trace Gases", fontsize=16,backgroundcolor = :transparent)
ax1 = Axis(f[1,1],aspect = DataAspect(),  title="Carbon-I Global Mode (200kg/hr)",xlabel = "meters",ylabel = "meters") 
sample=convert(Int,330/5)
hm = CairoMakie.heatmap!(ax1, dim_x[30:sample:end],dim_y[30:sample:end],  scene_global[30:sample:end,30:sample:end]./noise_global, colormap=:OrRd_9)
#hideydecorations!(ax1)
#hidexdecorations!(ax1)
CairoMakie.xlims!(-200,700)
CairoMakie.ylims!(-200,700)

cb = Colorbar(f[1,2], hm,  label = L"\text{\Delta \Omega(CH_4)/1\sigma}",  labelpadding = 10, ticklabelpad = 2, width = 10)
f
save("plots/PixelEnhancements300M.pdf",f)

### All in One:
function plotAll()
    f = Figure(resolution=(450,1000), title="Vertical Optical Depth of Trace Gases", fontsize=16,backgroundcolor = :transparent)
    ax1 = Axis(f[1,1],aspect = DataAspect(),  title="LES 5m simulation of a 200kg/hr point source") 
    sample=convert(Int,330/5)
    hm = CairoMakie.heatmap!(ax1, dim_x,dim_y,  C_sum, colormap=:OrRd_9)
    hideydecorations!(ax1)
    hidexdecorations!(ax1)
    CairoMakie.xlims!(-200,700)
    CairoMakie.ylims!(-200,700)

    cb = Colorbar(f[1,2], hm,  label = L"\text{\Delta \Omega(CH_4)} (mol/m²)",  labelpadding = 10, ticklabelpad = 2, width = 10,height = Relative(0.85))
    #save("plots/PixelLES.pdf",f)

    ax2 = Axis(f[2,1], aspect = DataAspect(),title="Carbon-I Target Mode (200kg/hr)") 
    sample=convert(Int,35/5)
    hm = CairoMakie.heatmap!(ax2, dim_x[1:sample:end],dim_y[1:sample:end],  scene_target[1:sample:end,1:sample:end]./noise_target, colormap=:OrRd_9)
    hideydecorations!(ax2)
    hidexdecorations!(ax2)
    #CairoMakie.ylims!(1e-5,1)
    CairoMakie.xlims!(-200,700)
    CairoMakie.ylims!(-200,700)
    cb2 = Colorbar(f[2,2], hm,  label = L"\text{\Delta \Omega(CH_4)/1\sigma}",  labelpadding = 10, ticklabelpad = 2, width = 10)
    #f
    #save("plots/PixelEnhancements30M.pdf",f)


    #f = Figure(resolution=(450,400), title="Vertical Optical Depth of Trace Gases", fontsize=16,backgroundcolor = :transparent)
    ax3 = Axis(f[3,1],aspect = DataAspect(),  title="Carbon-I Global Mode (200kg/hr)") 
    sample=convert(Int,330/5)
    start = 1
    hm = CairoMakie.heatmap!(ax3, dim_x[start:sample:end],dim_y[start:sample:end],  scene_global[start:sample:end,start:sample:end]./noise_global, colormap=:OrRd_9)
    hideydecorations!(ax3)
    hidexdecorations!(ax3)
    CairoMakie.xlims!(-200,700)
    CairoMakie.ylims!(-200,700)

    cb3 = Colorbar(f[3,2], hm,  label = L"\text{\Delta \Omega(CH_4)/1\sigma}",  labelpadding = 10, ticklabelpad = 2, width = 10)
    f
end
f = with_theme(plotAll, theme_black())
#CairoMakie.save("../../plots/ch4_fit_citb_dark.pdf", f)
save("plots/PixelEnhancementsAll_dark.pdf",f)


