function [Projs,param] = Cal_Bproj_2type(para, xdata, ydata)

Z_arr	   = xdata.Z_arr;
Res        = xdata.Res;
halfWidth  = xdata.halfWidth;
% step_sz    = xdata.step_sz;

model      = xdata.model;
angles     = xdata.angles;
atom       = xdata.atoms;
num_atom   = numel(atom);
atom_type_num = numel(unique(atom));

[N1,N2,~] = size(ydata);
num_pj=size(angles,1);
% N_s = 2*halfWidth+1;


fixedfa = reshape( make_fixedfa_man([N1,N2], Res, Z_arr), [N1,N2] );
% max_fa = max(abs(fixedfa(:)));
model = model/Res;

dtype = 'single';
X_rot = zeros(num_pj,num_atom,dtype);
Y_rot = zeros(num_pj,num_atom,dtype);
Z_rot = zeros(num_pj,num_atom,dtype);

for i=1:num_pj
  R1 = MatrixQuaternionRot([0 0 1],angles(i,1));  
  R2 = MatrixQuaternionRot([0 1 0],angles(i,2));
  R3 = MatrixQuaternionRot([1 0 0],angles(i,3));  
  R   = (R1*R2*R3)';
      
  rotCoords = R*model;
  X_rot(i,:) = rotCoords(1,:);
  Y_rot(i,:) = rotCoords(2,:);
  Z_rot(i,:) = rotCoords(3,:);
end

[X_crop,Y_crop] = ndgrid( -halfWidth:halfWidth, -halfWidth:halfWidth);      
Z_crop = -halfWidth:halfWidth;


%Projs   = zeros(N1,N1,num_pj);
h       = para(1,:)/para(1,1);
b       = (pi*Res)^2 ./ para(2,:);

num_atom_type = zeros(atom_type_num,1);
for j=1:atom_type_num
    num_atom_type(j) = nnz(atom==j);
end

% t = (step_sz/max_fa^2/num_pj/N1^2) ./ num_atom_type;
% sum_pj = sqrt(sum(ydata(:).^2));

Grad    = zeros(N1,N2,num_pj, 3);
% grad_h  = zeros(N1,N2,num_pj, 3);
% grad_b  = zeros(N1,N2,num_pj, 3);
% tic
for i=1:num_pj

    for j=1:atom_type_num
        atom_type_j = atom==j;
        
        X_cen = reshape(X_rot(i,atom_type_j), [1,1,num_atom_type(j)]);
        Y_cen = reshape(Y_rot(i,atom_type_j), [1,1,num_atom_type(j)]);
        Z_cen = reshape(Z_rot(i,atom_type_j), [1,1,num_atom_type(j)]);
        
        X_round = round(X_cen);
        Y_round = round(Y_cen);
        Z_round = round(Z_cen);
        
        l2_xy = bsxfun(@plus, X_crop, X_round - X_cen).^2 + ...
                bsxfun(@plus, Y_crop, Y_round - Y_cen).^2;
        
        l2_z  = bsxfun(@plus, Z_crop, Z_round - Z_cen).^2    ;        
        
        pj_j   = bsxfun(@times, exp(-l2_xy* b(j) ), sum(exp(-l2_z* b(j))) );
        pj_j_h = h(j)*pj_j;
%         bj_j   = pj_j_h.*l2_xy  + h(j)* bsxfun(@times, exp(-l2_xy* b(j) ), sum(l2_z.*exp(-l2_z* b(j))) );
        
        
        for k=1:num_atom_type(j)
            indx = X_round(k) + (-halfWidth:halfWidth) + round((N1+1)/2);
            indy = Y_round(k) + (-halfWidth:halfWidth) + round((N2+1)/2);

            Grad(indx,indy,i,j)   = Grad(indx,indy,i,j)   + pj_j_h(:,:,k);
%             grad_h(indx,indy,i,j) = grad_h(indx,indy,i,j) + pj_j(:,:,k);
%             grad_b(indx,indy,i,j) = grad_b(indx,indy,i,j) + bj_j(:,:,k);
        end
    end
end

Projs = sum(Grad,4);
for i=1:num_pj
    Projs(:,:,i) = (( my_ifft( my_fft(Projs(:,:,i)) .* (fixedfa) ) ));
end
% res = (Projs - ydata);
% fprintf('%d.f = %.5f\n', iter, norm( res(:) ) /sum_pj );


k=sum(Projs(:).*ydata(:))/sum(Projs(:).^2);
% k=sum(ydata(:))/sum(Projs(:));
Projs=Projs*k;


% param=[h(:)'; (pi*Res)^2 ./ b(:)'];
param=[k*h; (pi*Res)^2 ./ b(:)'];
end