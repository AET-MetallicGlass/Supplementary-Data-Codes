function [Projs,params] = gradient_XYZ_man(para, xdata, ydata)
fprintf('\nHB gradient algorithm\n');

Z_arr	   = xdata.Z_arr;
Res        = xdata.Res;
halfWidth  = xdata.halfWidth;
iterations = xdata.iterations;

model      = xdata.model;
model_ori  = xdata.model_ori;
angles     = xdata.angles;
atom       = xdata.atoms;
num_atom   = numel(atom);
step_sz    = xdata.step_sz;

if numel(step_sz) == 1
    step_sz = step_sz.*ones(1,iterations);
elseif numel(step_sz) < iterations
    step_sz(numel(step_sz)+1:iterations) = step_sz(numel(step_sz));
end

[N1,~,num_pj] = size(ydata);

% atom_type1 = xdata.atom_type1;
% atom_type2 = xdata.atom_type2;
% atom_type3 = xdata.atom_type3;
% X_rot = xdata.X_rot;
% Y_rot = xdata.Y_rot;
% Z_rot = xdata.Z_rot;

fixedfa = reshape( make_fixedfa_man(N1, Res, Z_arr), [N1,N1] );
model = model/Res;
model_ori = model_ori/Res;
%
for i=1:num_pj
    ydata(:,:,i) = max(real( my_ifft( my_fft( ydata(:,:,i) ) ./fixedfa )),0);
end

% l2_type1 =  xdata.l2_type1;
% l2_type2 =  xdata.l2_type2;
% l2_type3 =  xdata.l2_type3;

% h1 = para(1,1);
% h2 = para(1,2);
% h3 = para(1,3);
%
% b1 = para(2,1);
% b2 = para(2,2);
% b3 = para(2,3);


[X_crop,Y_crop] = ndgrid( -halfWidth:halfWidth, -halfWidth:halfWidth);
Z_crop = -halfWidth:halfWidth;

[xx,yy,zz] = ndgrid( -halfWidth:halfWidth, -halfWidth:halfWidth, -halfWidth:halfWidth);
R2 = xx.^2 + yy.^2 + zz.^2;
sum_R2 = sum(R2(:));
clear xx yy zz R2;

%Y_crop = reshape(Y_crop, [1,2*halfWidth+1,2*halfWidth+1]);
%X_crop = reshape(X_crop, [1,2*halfWidth+1,2*halfWidth+1]);
%Z_crop = reshape(Z_crop, [1,2*halfWidth+1]);


h = zeros(1,1,num_atom);
b = zeros(1,1,num_atom);
for k=1:size(para,2)
    h(atom==k) = para(1,k);
    b(atom==k) = para(2,k);
end
% if size(para,2)==3
%     for k=1:3
%         h(atom==k) = para(1,k);
%         b(atom==k) = para(2,k);
%     end
% elseif size(para,2)==num_atom
%     h(:) = para(1,:);
%     b(:) = para(2,:);
% else
%     disp('error')
%     return
% end
b = (Res*pi)^2./b;

%num_atom_type = [nnz(atom==1),nnz(atom==2),nnz(atom==3)];
sum_pj = sqrt(sum(ydata(:).^2));

N_s = 2*halfWidth+1;
index = zeros(2,num_pj, num_atom);

grad_x_set = zeros(N_s,N_s,num_pj, num_atom,'single');
grad_y_set = zeros(N_s,N_s,num_pj, num_atom,'single');
grad_z_set = zeros(N_s,N_s,num_pj, num_atom,'single');

% for i=1:num_pj
%     X_round = round(X_rot(i,:));
%     Y_round = round(Y_rot(i,:));
%     
%     for k=1:num_atom
%         indx = X_round(k) + (-halfWidth:halfWidth) + N1/2 +1;
%         indy = Y_round(k) + (-halfWidth:halfWidth) + N1/2 +1;
%         [yy,xx] = meshgrid(indy, indx);
%         index_set( (i-1)*N_s^2 + 1 : i*N_s^2, k) = sub2ind( [N1,N1,num_pj], xx(:),yy(:),i*ones(N_s^2,1) );
%         %temp = false(N1,N1);temp(indx,indy)=1;
%         %logical_set((i-1)*N1^2 + 1 : i*N1^2, k) = logical(temp(:));
%         
%     end
% end

h = (single(h));
b = (single(b));

%index_set = (index_set);
X_crop = (single(X_crop));
Y_crop = (single(Y_crop));
Z_crop = (single(Z_crop));

X_ori = reshape( model_ori(1,:), [1,1,num_atom] );
Y_ori = reshape( model_ori(2,:), [1,1,num_atom] );
Z_ori = reshape( model_ori(3,:), [1,1,num_atom] );
X     = reshape( model(1,:), [1,1,num_atom] );
Y     = reshape( model(2,:), [1,1,num_atom] );
Z     = reshape( model(3,:), [1,1,num_atom] );

if isfield(xdata, 'saveOutput'), radius = xdata.radius;
else  radius = 1;
end
scale=radius/Res;

model = zeros(size(model,1),size(model,2),numel(step_sz));

for iter=1:iterations
    dt = (num_pj*step_sz(iter) / (8*sum_R2*num_pj*N1) ) ./ ( h.^2 .* b.^2 );
    
    Projs = zeros(N1,N1,num_pj,'single');

    %tic
    for i=1:num_pj
        
        RM1 = MatrixQuaternionRot([0 0 1],angles(1,i));
        RM2 = MatrixQuaternionRot([0 1 0],angles(2,i));
        RM3 = MatrixQuaternionRot([1 0 0],angles(3,i));
        R   = RM1*RM2*RM3;
        model_rot = R'*[X(:)';Y(:)';Z(:)'];
        
        X_cen = reshape(model_rot(1,:), [1,1,num_atom]);
        Y_cen = reshape(model_rot(2,:), [1,1,num_atom]);
        Z_cen = reshape(model_rot(3,:), [1,1,num_atom]);
        
        X_round = round(X_cen);
        Y_round = round(Y_cen);
        Z_round = round(Z_cen);
        
        Dx = bsxfun(@plus, X_crop, X_round - X_cen);
        Dy = bsxfun(@plus, Y_crop, Y_round - Y_cen);
        Dz = bsxfun(@plus, Z_crop, Z_round - Z_cen);        
        
        l2_xy = Dx.^2 + Dy.^2;
        l2_z  = Dz.^2;

        l2_xy_b     = bsxfun(@times, l2_xy, b);
        l2_z_b      = bsxfun(@times, l2_z , b);
        exp_l2_z_b  = exp(-l2_z_b);
        exp_l2_xy_b = exp(-l2_xy_b);
        
        pj_j     = bsxfun(@times, exp_l2_xy_b, sum(exp_l2_z_b) );
        pj_j_h   = bsxfun(@times, pj_j , h);
        pj_j_b   = bsxfun(@times, pj_j_h , b);
                
        R2_Dx = ( R(1,1)*Dx + R(1,2)*Dy ) .* pj_j_b;
        R2_Dy = ( R(2,1)*Dx + R(2,2)*Dy ) .* pj_j_b;
        R2_Dz = ( R(3,1)*Dx + R(3,2)*Dy ) .* pj_j_b;
        
        Dz_exp_l2_z_b = Dz .* exp_l2_z_b;
        sum_Dz_exp    = bsxfun(@times, sum( Dz_exp_l2_z_b ), exp_l2_xy_b);
        
        sum_Dz_hb     = bsxfun(@times, sum_Dz_exp, h.*b);        
        
        xj_j = R2_Dx + R(1,3) * sum_Dz_hb;
        yj_j = R2_Dy + R(2,3) * sum_Dz_hb;
        zj_j = R2_Dz + R(3,3) * sum_Dz_hb;
                
        for k=1:num_atom
            indx = X_round(k) + (-halfWidth:halfWidth) + N1/2 +1;
            indy = Y_round(k) + (-halfWidth:halfWidth) + N1/2 +1;

            Projs(indx,indy,i) = Projs(indx,indy,i) + pj_j_h(:,:,k);
            %grad_h(:,:,i,k) = pj_j(:,:,k);            
            index(:,i,k) = [X_round(k);Y_round(k)];            
            
            grad_x_set(:,:,i,k) = xj_j(:,:,k);
            grad_y_set(:,:,i,k) = yj_j(:,:,k);
            grad_z_set(:,:,i,k) = zj_j(:,:,k);
        end
        
        %Projs(:,:,i) = real( my_ifft( my_fft( Projs(:,:,i)) .* fixedfa ) );
    end
    
    
    % for i=1:num_pj
    %     Projs(:,:,i) = real(( my_ifft( my_fft(Projs(:,:,i)) .* fixedfa ) ));
    % end
    res = Projs - ydata;
    fprintf('%d.f = %.5f\n', iter, sqrt(sum(res(:).^2)) /sum_pj );
    grad_x = zeros(1,1,num_atom,'single'); 
    grad_y = zeros(1,1,num_atom,'single'); 
    grad_z = zeros(1,1,num_atom,'single'); 
    for k=1:num_atom
        for i=1:num_pj
            indx = index(1,i,k) + (-halfWidth:halfWidth) + N1/2 +1;
            indy = index(2,i,k) + (-halfWidth:halfWidth) + N1/2 +1;
            grad_x(k) = grad_x(k) + sum(sum( res(indx,indy,i) .* grad_x_set(:,:,i,k) ));
            grad_y(k) = grad_y(k) + sum(sum( res(indx,indy,i) .* grad_y_set(:,:,i,k) ));
            grad_z(k) = grad_z(k) + sum(sum( res(indx,indy,i) .* grad_z_set(:,:,i,k) ));                        
        end
    end
    %if mod(iter,400)==0, dt = dt/2; end
%     fprintf('[%f,%f,%f\n',[max(abs(dt(:).*grad_x(:))),max(abs(dt(:).*grad_y(:))),max(abs(dt(:).*grad_z(:)))] );
    dx = dt.*grad_x; dx=min(dx,.5); dx = max(-.5, dx);
    dy = dt.*grad_y; dy=min(dy,.5); dy = max(-.5, dy);
    dz = dt.*grad_z; dz=min(dz,.5); dz = max(-.5, dz);
    X = X - dx;
    Y = Y - dy;
    Z = Z - dz;
        
    diff_X = X-X_ori;
    diff_Y = Y-Y_ori;
    diff_Z = Z-Z_ori;
    diff_norm = sqrt( diff_X.^2 + diff_Y.^2 + diff_Z.^2 );
    index_3 = diff_norm>scale;
    X(index_3) = X_ori(index_3) + scale*diff_X(index_3) ./ diff_norm(index_3);
    Y(index_3) = Y_ori(index_3) + scale*diff_Y(index_3) ./ diff_norm(index_3);
    Z(index_3) = Z_ori(index_3) + scale*diff_Z(index_3) ./ diff_norm(index_3);
    
    model(:,:,iter) = [X(:),Y(:),Z(:)]'*Res;
    %toc
end
%%
for i=1:num_pj
    Projs(:,:,i) = max( 0, real( my_ifft( my_fft(Projs(:,:,i)) .* fixedfa ) ) );
end

params = model;
end