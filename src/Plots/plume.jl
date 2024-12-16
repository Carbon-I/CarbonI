les = Dataset("/net/fluo/data2/groupMembers/siraput/Data_NewSatRetrieval/LES_P_T_W/C_conc_30817.h5")

C_ = les["C"][:]
dims = size(C_)
area = (40,40)
C_3D = les["C"];
sub = C_3D[:,100:750,:]
C2D = sum(sub, dims=1)[1,:,:];
dims = size(C2D)
area = (12,12)
area = (44,60)
C_sub = zeros(dims[1]-area[1], dims[2]-area[2])
for i=1:area[1], j=1:area[2]
    #@show size(C_[100, i:end-area[1]+i-1,i:end-area[2]+i-1]), size(C_sub)
    C_sub += (C2D[i:end-area[1]+i-1,j:end-area[2]+j-1])
end
C_sub /= area[1]*area[2]
CC = C_sub[1:6:end,1:6:end]';

CC[CC.<0.5e17].=NaN

bck = la[3000:3250,2000:2250,5]
cloud = deepcopy(bck);
cloud[bck.<0.07].=NaN
bck[bck.>0.07].=NaN
la[3000:3250,2000:2250,5]