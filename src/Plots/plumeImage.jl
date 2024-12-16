using NCDatasets
using StatsBase

les = Dataset("/net/fluo/data2/groupMembers/siraput/Data_NewSatRetrieval/LES_P_T_W/C_conc_30817.h5")

C3D = les["C"][:];
C_sum = sum(C3D, dims=1)[1,:,:];


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

C1 = C_sum'[:,100:400];
C6 = downsample(C_sum'[:,100:400],6);
C60 = downsample(C_sum'[:,100:400],60);
x1 = range(0, size(C1,1)*5, size(C1,1))
y1 = range(0, size(C1,2)*5, size(C1,2))
 
x6 = range(0, (size(C6,1)-1)*30, size(C6,1))
y6 = range(0, (size(C6,2)-1)*30, size(C6,2))
x60 = range(0, (size(C60,1)-1)*300, size(C60,1))
y60 = range(0, (size(C60,2)-1)*300, size(C60,2))

p1 = Plots.heatmap(y1,x1,C1,xaxis=nothing, yaxis=nothing, cmap=:gist_yarg, aspect_ratio=:equal, clim=(0,0.5e18), title="5m resolution", colorbar=:false); xlims!(0,1190); ylims!(0,1005)
p2 = Plots.heatmap(y6,x6,C6,xaxis=nothing, yaxis=nothing,cmap=:gist_yarg, aspect_ratio=:equal, clim=(0,0.5e18), title="30m resolution", colorbar=:false); xlims!(0,1190); ylims!(0,1005)
p3 = Plots.heatmap(y60,x60,C60,xaxis=nothing, yaxis=nothing,cmap=:gist_yarg, aspect_ratio=:equal, clim=(0,0.5e18), title="300m resolution",colorbar=:false); xlims!(0,1190); ylims!(0,1005)

l = @layout [a b c]

Plots.plot(p1, p2, p3, layout = l)
Plots.plot!(size=(950,300))
savefig("/home/cfranken/PixelEnhancements_map.pdf")
savefig("/home/cfranken/PixelEnhancements_map.png")

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