function new_atomtype = My_reClass_3atoms(DataMatrix,curr_model,curr_atomtype,Rad,halfSize,SPHyn)

[xX,yY,zZ] = ndgrid(-halfSize:halfSize,-halfSize:halfSize,-halfSize:halfSize);
SphereInd = find(xX.^2+yY.^2+zZ.^2 <=(halfSize+0.5)^2);

if SPHyn==1
    useInd = SphereInd;
else
    useInd = 1:length(xX(:));
end

new_atomtype = curr_atomtype;

for i=1:size(curr_model,2)
    curr_atompos = curr_model(:,i);
    Dist = sqrt(sum((curr_model - repmat(curr_atompos,[1 size(curr_model,2)])).^2,1));
    
    BallInd = find(Dist~=0 & Dist < Rad);
    
    BallType = curr_atomtype(BallInd);
    BallModel = curr_model(:,BallInd);
    
    [avatom1]= compute_average_atom_from_vol(DataMatrix,BallModel,find(BallType==1),halfSize);
    [avatom2]= compute_average_atom_from_vol(DataMatrix,BallModel,find(BallType==2),halfSize);
    [avatom3]= compute_average_atom_from_vol(DataMatrix,BallModel,find(BallType==3),halfSize);
    
    [new_atomtype(i)] = get_atomtype_vol_threeatom_noOutPut_useInd(DataMatrix, curr_atompos, avatom1, avatom2, avatom3, halfSize,useInd);
    
end

end