clear all; close all
%%
addpath('data\')
addpath('src\')
%% read projections after BM3D and tilt angles
projections = importdata('CoPdPt_Proj_Dose1p7e4_BM3D.mat');
angles      = importdata('CoPdPt_Angle_Dose1p7e4.mat');
%% input
rotation       = 'ZYX';  % Euler angles setting ZYZ
dtype          = 'single';
projections_refined = cast(projections,dtype);
angles_refined      = cast(angles,dtype);

% compute normal vector of rotation matrix
matR = zeros(3,3);
if length(rotation)~=3
    disp('rotation not recognized. Set rotation = ZYX\n'); rotation = 'ZYX';
end
for i=1:3
    switch rotation(i)
        case 'X',   matR(:,i) = [1;0;0];
        case 'Y',   matR(:,i) = [0;1;0];
        case 'Z',   matR(:,i) = [0;0;1];
        otherwise,  matR = [0,0,1;
                0,1,0;
                1,0,0];
            disp('Rotation not recognized. Set rotation = ZYX');
            break
    end
end
vec1 = matR(:,1); vec2 = matR(:,2); vec3 = matR(:,3);

% extract size of projections & num of projections
[dimx, dimy, Num_pj] = size(projections_refined);
%% parameter
step_size      = 2;  
iterations     = 300;
dimz           = dimx;
positivity     = true;
%% rotation matrix
Rs = zeros(3,3,Num_pj, dtype);
for k = 1:Num_pj
    phi   = angles_refined(k,1);
    theta = angles_refined(k,2);
    psi   = angles_refined(k,3);
    
    % compute rotation matrix R w.r.t euler angles {phi,theta,psi}
    rotmat1 = MatrixQuaternionRot(vec1,phi);
    rotmat2 = MatrixQuaternionRot(vec2,theta);
    rotmat3 = MatrixQuaternionRot(vec3,psi);
    R =  single(rotmat1*rotmat2*rotmat3)';
    Rs(:,:,k) = R;
end

%% Run RESIRE with multi-GPU
fprintf('\nResire code:(multi-GPU)\n');
dim_ext = [dimx,dimy,dimz];
tic
[rec] = RT3_film_multiGPU( (projections_refined), (Rs), dim_ext, ...
    iterations, (step_size) , (positivity) );

%% Cal. Projections
ref_projs = calculate3Dprojection_multiGPU(single(rec), Rs);
figure();img(ref_projs,'colormap','gray')

%% Output
OBJ.reconstruction=rec;
% save("Resire_CoPdPt_Dose1p7e4_rec.mat","OBJ",'-v7.3')
% save("Resire_CoPdPt_Dose1p7e4_CalProj.mat","ref_projs",'-v7.3')
%% Visualization of rec
Mask=importdata('Resire_CoPdPt_Mask.mat');

for i=1:346
    Dsetvol(:,:,i)=rec(:,:,i).*Mask(:,:,i);
end
figure
img(Dsetvol(:,:,40:306),[],'colormap','gray')
return
%% Visualization of Projs.
load('Resire_CoPdPt_Dose1p7e4_CalProj.mat')
load('Proj_Dose1p7e4_PoissonNoise_BeforeBM3D.mat')
%
projections_BM3D=projections;
figure();
img(Proj_Dose1p7e4_PoissonNoise(30:end-30,30:end-30,:),[],projections_BM3D(30:end-30,30:end-30,:),[],ref_projs(30:end-30,30:end-30,:),[],'colormap','gray')
%%