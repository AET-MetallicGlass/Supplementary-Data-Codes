%% Main_classification %%
% perform both global and local k-mean on the summed intensity of 9*9*9
% volume around traced atom position

addpath('src/')
addpath('input/')

% add the path to load reconstruction volume, you can comment it and move
% the reconstruction into input folder
addpath('../3_Final_reconstruction_volume/') 

% read in files: traced atomic position and reconstruction volume
% due to the manually chceck and iteration process between the position 
% refinement and classification, here only shows the last round results 
% which could be some deviation from the results achieved by 
% Main_polynomial_tracing.m
new_model   = importdata('traced_model_inPixel.mat');
Dsetvol     = importdata('MG_reconstruction_volume.mat');

%% upsampling the reconstruction matrix by 3*3*3 by linear interpolation
% better to run at super conputer since the size of the interpolated
% volume will be larger than 16G

xx = (1:size(Dsetvol,1)) - round((size(Dsetvol,1)+1)/2);
yy = (1:size(Dsetvol,2)) - round((size(Dsetvol,2)+1)/2);
zz = (1:size(Dsetvol,3)) - round((size(Dsetvol,3)+1)/2);

xxi = ((3*xx(1)):(xx(end)*3))/3;
yyi = ((3*yy(1)):(yy(end)*3))/3;
zzi = ((3*zz(1)):(zz(end)*3))/3;

xxi = xxi(3:end); yyi = yyi(3:end); zzi = zzi(3:end);

[Y,X,Z] = meshgrid(yy,xx,zz);
[Yi,Xi,Zi] = meshgrid(yyi,xxi,zzi);

Dsetvol = interp3(Y,X,Z,Dsetvol,Yi,Xi,Zi,'spline',0);

clear Xi Yi Zi
FinalVol = My_paddzero(Dsetvol,size(Dsetvol)+20);

FinalVol_single = single(FinalVol);
clear FinalVol Dsetvol
%% apply global k-mean classification on the reconstruction with traced atom position

% set parameters
% the half size of volume around atom to calculate intensity, halfsize = 3
% means to sum a 7*7*7 box around the atom to get the intensity
classify_info = struct('Num_species', 3,  'halfSize',  3,  'plothalfSize',  1, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  70);

% calculate the atomic position with upsampled volume
new_model_L = (new_model +2).*3;
% apply k-mean classification on the reconstruction
[atom_model, global_class_atomtype] = initial_class_kmean_sub(...
    FinalVol_single, new_model_L, classify_info);
% apply function plot_class_hist() to achieve the histogram information peak_info
% please see the descriptions in subfunction to get more details
[peak_info_global_classfication,~] = plot_class_hist(...
    FinalVol_single, atom_model, global_class_atomtype, classify_info);
%% apply local k-mean classification on the reconstruction by the results of global k-mean

temp_class_atomtype = global_class_atomtype;

classify_info.Radius = 10/0.347*3; % radius is 10A

% when there are 5 iterations with same number of atoms flipped (back and
% forth), the interation will be stopped
[atom_model, local_class_atomtype] = local_class_kmean_sub(...
    FinalVol_single, atom_model, temp_class_atomtype, classify_info);

% apply function plot_class_hist() to achieve the histogram information peak_info
% please see the descriptions in subfunction to get more details
[peak_info_local_classfication,~] = plot_class_hist(...
    FinalVol_single, atom_model, local_class_atomtype, classify_info);

save('output/localC_res.mat','local_atomtype');
