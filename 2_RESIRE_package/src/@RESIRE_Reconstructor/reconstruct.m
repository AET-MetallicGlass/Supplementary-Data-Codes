%% RESIRE reconstruction

%%inputs:
%%  InputProjections - measured projections
%%  Support - region of 1's and 0's defining where the reconstruction can exist 
%%    (voxels that are 0 in support will be set to 0 in reconstruction)
%%  sum_rot_pjs - the sum of constant terms for the gradient descent: sum(A^T(b)),
%%      A is the back projection operator, b is measured projection
%%  xj, yj, zj - different dimension for rotated coordinates with oversampling
%%  Rot_x, Rot_y - x and y dimension for rotated coordinates without oversampling
%%  step_size - the step size for gradient descent
%%  dtype - the type of output results, set as 'single' to save the memory
%%  monitor_R - the flag for recording real space error

%%outputs:
%%  reconstruction - reconstruction after iteration
%%  errR - Real space error (averaged L1 norm between the measured projections and backprojections)
%%  Rarr_record - recorded L1 norm error array for all projections
%%  Rarr2_record - recorded L2 norm error array for all projections

function obj = reconstruct(obj)
projections = obj.InputProjections;
Num_pj      = obj.NumProjs;
step_size   = obj.step_size;
dimx        = obj.Dim1;
dimy        = obj.Dim2;
dtype       = obj.dtype;
Rot_x       = obj.Rot_x;
Rot_y       = obj.Rot_y;
iterations  = obj.numIterations;

xj = double(obj.xj); yj = double(obj.yj); zj = double(obj.zj);

sum_rot_pjs = obj.sum_rot_pjs; 

rec = zeros(dimy,dimx,dimy,dtype);
rec_big = zeros(obj.n2_oversampled,obj.n1_oversampled,obj.n2_oversampled,dtype);
ind_V	= My_volumn_index(size(rec_big),size(rec));
if obj.initial_model == 1
    rec = obj.Support;
    rec_big(ind_V(1,1):ind_V(1,2), ind_V(2,1):ind_V(2,2), ind_V(3,1):ind_V(3,2) ) = rec;
end
dt      = (step_size/Num_pj/dimx);

fprintf('RESIRE: Reconstructing... \n\n');

if obj.monitor_R == 1
    monitorR_loopLength = obj.monitorR_loopLength;
    errR_arr = zeros(1,floor(iterations./monitorR_loopLength));
    Rarr_record = zeros(Num_pj,floor(iterations./monitorR_loopLength));
    Rarr2_record = zeros(Num_pj,floor(iterations./monitorR_loopLength));
end

% flag for using parallel calculation
if obj.use_parallel
  parforArg = Inf;
else
  parforArg = 0;
end

for iter=1:iterations
    recK  = double(my_fft(rec_big));    
    
    % compute rotated projections via Fourier Slice Theorem
    pj_cal = splinterp3(recK, yj, xj, zj);
    pj_cal = real(fftshift(ifft2(ifftshift(pj_cal))));
    pj_cal = croppedOut(pj_cal, [dimx,dimy,Num_pj] );
    
    % compute R factor
    if obj.monitor_R == 1 && mod(iter,monitorR_loopLength) == 0
        Rarr = zeros(Num_pj,1);
        Rarr2 = zeros(Num_pj,1);
        for i=1:Num_pj
            pj = projections(:,:,i); 
            proj_i = pj_cal(:,:,i);
            Rarr(i) = sum(sum( abs(proj_i - pj) ))/ sum(abs(pj(:)));
            Rarr2(i) = norm( proj_i - pj, 'fro' )/ norm( pj, 'fro' );
        end
        errR  = mean(Rarr);
        errR2 = mean(Rarr2);
        fprintf('RESIRE: Iteration %d. Rfactor=%.4f, R2factor=%.4f \n',iter, errR, errR2);
        errR_arr(iter./monitorR_loopLength) = errR;
        Rarr_record(:,iter./monitorR_loopLength) = Rarr;
        Rarr2_record(:,iter./monitorR_loopLength) = Rarr2;
    else
        fprintf('RESIRE: Iteration %d \n',iter);
    end    
    
    % compute gradient & apply gradient descent
    grad = -sum_rot_pjs;
    parfor (k = 1:Num_pj,parforArg)
%     for k = 1:Num_pj
        rot_pj_cal  = splinterp2(pj_cal(:,:,k), double(Rot_y(:,:,:,k)), double(Rot_x(:,:,:,k)));
        grad = grad + rot_pj_cal;
    end
    rec = rec - dt*grad;
    rec = max(0,rec);
    
    % flag for saving temporary reconstruction
    if obj.save_temp == 1 && mod(iter,obj.save_loopLength) == 0
        temp_filename = [obj.saveFilename,num2str(iter),'.mat'];
        save(temp_filename, 'rec');
    end
    
    rec_big( ind_V(1,1):ind_V(1,2), ind_V(2,1):ind_V(2,2), ind_V(3,1):ind_V(3,2) ) = rec;
end
if obj.monitor_R == 1
    obj.errR            = errR;
    obj.Rarr_record     = Rarr_record;
    obj.Rarr2_record    = Rarr2_record;
end
obj.reconstruction = rec;
end