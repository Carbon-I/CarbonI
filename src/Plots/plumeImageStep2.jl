using NCDatasets
using StatsBase, LaTeXStrings, CairoMakie

les = Dataset("/net/fluo/data1/ftp/XYZT_LES_/Data_Train_Val_Test/original_train_set.h5")
#les = Dataset("/net/fluo/data2/groupMembers/siraput/Data_NewSatRetrieval/LES_P_T_W/C_conc_27777.h5")
#les = Dataset("/net/fluo/data2/groupMembers/siraput/Data_NewSatRetrieval/LES_P_T_W/run10_C/run10_C_26902.h5")
C3D = les.group["train_set"]["C_2D"][:,:,10,2];
C_sum = sum(C3D, dims=1)[1,:,:];

# From Siraput, unit conversion into ppm-m (in /net/fluo/data1/ftp/XYZT_LES_/code/helper_functions.py)
c2ppm = (500.0/3.57e-4)
c2mol = c2ppm*44/1e6
C_sum = C3D*10*c2mol # Use 200kg/hr (ref is about 5g/s
function downsample(C_sum, factor)    
    C_sum_ds = zeros(size(C_sum,1)÷factor, size(C_sum,2)÷factor);
    for i=1:size(C_sum,1)÷factor
        for j=1:size(C_sum,2)÷factor
            C_sum_ds[i,j] = mean(C_sum[(i-1)*factor+1:i*factor, (j-1)*factor+1:j*factor]);
        end
    end 
    return C_sum_ds
end

PlumeStats = zeros(100, 3);
for factor = 1:100
    shift = max(1,factor ÷ 2);
    
    CC  = downsample(C_sum',factor);
    CC1 = downsample(C_sum[shift:end,:]',factor);
    CC2 = downsample(C_sum[:,shift:end]',factor);
    CC3 = downsample(C_sum[shift:end,shift:end]',factor);
    
    pix1 = [sort(A)[end] for A in (CC[:],CC1[:],CC2[:],CC3[:])];
    #@show pix1
    pix3 = [sort(A)[end-2] for A in (CC[:],CC1[:],CC2[:],CC3[:])];
    pix10 = [sort(A)[end-9] for A in (CC[:],CC1[:],CC2[:],CC3[:])];
    #PlumeStats[factor,1] = maximum(CC);
    PlumeStats[factor,1] = mean(pix1);
    PlumeStats[factor,2] = mean(pix3);
    PlumeStats[factor,3] = mean(pix10);
end



C1 = C_sum';
C6 = downsample(C_sum',6);
C60 = downsample(C_sum',60);
x1 = range(0, size(C1,1)*5, size(C1,1))
y1 = range(0, size(C1,2)*5, size(C1,2))
 
x6 = range(0, (size(C6,1)-1)*30, size(C6,1))
y6 = range(0, (size(C6,2)-1)*30, size(C6,2))
x60 = range(0, (size(C60,1)-1)*300, size(C60,1))
y60 = range(0, (size(C60,2)-1)*300, size(C60,2))



refColumn  = 4e19/6e23*100^2
noise_target  = 4 * sqrt(10) / 1900 * refColumn
noise_real_target = randn(size(C6')) * noise_target
f = Figure(resolution=(450,400), title="Vertical Optical Depth of Trace Gases", fontsize=16)
ax1 = Axis(f[1,1],  aspect = DataAspect(),title="200kg/hr at 30m GSD") 
hm = CairoMakie.heatmap!(ax1, C6'./noise_target, colormap=:OrRd_9)
hideydecorations!(ax1)
hidexdecorations!(ax1)
#CairoMakie.ylims!(1e-5,1)

cb = Colorbar(f[1,2], hm,  label = L"\text{\Delta \Omega(CH_4)/1\sigma}",  labelpadding = 10, ticklabelpad = 2, width = 10)
f
save("plots/PixelEnhancements30M.pdf",f)
noise_real_global = randn(size(C60')) * noise_target/sqrt(10)
f = Figure(resolution=(450,400), title="Vertical Optical Depth of Trace Gases", fontsize=16)
ax1 = Axis(f[1,1],aspect = DataAspect(),  title="200kg/hr at 300m GSD") 
hm = CairoMakie.heatmap!(ax1, C60'./(noise_target/sqrt(10)), colormap=:OrRd_9, aspect_ratio=:equal)
hideydecorations!(ax1)
hidexdecorations!(ax1)
#CairoMakie.xlims!(-100,800)
#CairoMakie.ylims!(1e-5,1)

cb = Colorbar(f[1,2], hm,  label = L"\text{\Delta \Omega(CH_4)/1\sigma}",  labelpadding = 10, ticklabelpad = 2, width = 10)
f
save("plots/PixelEnhancements300M.pdf",f)
#plot(5*(1:100), PlumeStats[:,1]/PlumeStats[60,1], yscale=:log10, label = "highest pixel enhancement")
plot(5*(1:100), PlumeStats[:,2]/PlumeStats[60,2], yscale=:log10, yticks=([0.5,1,2,4,7,10,15,20,30],[0.5,1,2,4,8,10,15,20,30]), label = "3rd highest pixel enhancement")
plot!(5*(1:100), PlumeStats[:,3]/PlumeStats[60,3], label = "10th highest pixel enhancement")
plot!([30,30],[0.5, 14], color=:black, linewidth=3, label="Carbon-I Target mode (CBE)") 
plot!([300,300],[0.5, 1], color=:gray, linewidth=3,label="Carbon-I Global mode (CBE)")
ylims!(0.5,100)
xlabel!("Footprint size (m)")
ylabel!("maximum enhancements (normalized to 300m)")
plot!(size=(600,450))
savefig("/home/cfranken/PixelEnhancements.pdf")