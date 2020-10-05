function [atomtype, Rs] = initial_class_L1norm(box_arr, mean_box, O_Ratio, halfSize, SPHyn )

atomtype = -1*ones(1,size(box_arr,4));
Rs = zeros(2,size(box_arr,4));

ds = 1/O_Ratio;

[XX,YY,ZZ] = ndgrid(-halfSize: ds :halfSize, ...
                    -halfSize: ds :halfSize, ...
                    -halfSize: ds :halfSize);
                
if SPHyn
    useInd = find( (XX.^2 + YY.^2 + ZZ.^2) <= (halfSize + 0.5*ds)^2 );
else
    useInd = 1:length(XX);
end

box2coordinates  = create_box(size(box_arr,1));

fit_param_init  = [0,  max(mean_box(:)),  0,  0,  0,  0.5,  0.5,    0.5,  0,  0,  0];
fixed           = [0,  0            ,  0,  0,  0,    0,    0,    0,    0,  0,    0];
lb              = [0,  0            , -2, -2, -2,    0,    0,    0, -pi,  0,  -pi];
ub              = [inf,inf          ,  2,  2,  2,  inf,  inf,  inf,  pi,  pi,  pi];

%perform the fit
if halfSize == 0
    fit_result = 0;
else
    [fit_result ] = fit_gauss3D_PD(fit_param_init, box2coordinates, mean_box, fixed, lb, ub);
end
fprintf('fitted constant = %4.02f\n',fit_result(1));

for ind = 1:size(box_arr,4)
    
    DataBox = box_arr(:,:,:,ind);
    
    ZeroBox = zeros(size(mean_box))+ fit_result(1);
%     ZeroBox = zeros(size(mean_box));
    %calculate squared deviation ( = non-normalized r-factor)
    R1 = sum(abs(DataBox(useInd) - ZeroBox(useInd)));
    R2 = sum(abs(DataBox(useInd) - mean_box(useInd)));
    Rs(1,ind) = R1;
    Rs(2,ind) = R2;
    
    if (R1 > R2)
        %this is an atom, add to density matrix
        atomtype(ind) = 1;
    else
        %not an atom, add background to density matrix
        atomtype(ind) = 0;
    end
        
end

fprintf('number of inserted atoms: \t %i atoms\n', sum(atomtype==1))
fprintf('number of  skipped atoms: \t %i atoms\n', sum(atomtype==0))
end