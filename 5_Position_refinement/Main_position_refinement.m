%% Main Position refinement 
% refinement of the atomic coordinates with existed model, type and
% reconstruction volume: minimize the error between the atomic coordinates 
% and the measured projections

addpath('src/')
addpath('input/')

% add the path to load measured projections and angles, you can comment it 
% and move the projections and angles into input folder
addpath('../1_Measured_data/') 

% read in files: measured projections and angles; 
% atomic position and types after classification
projections = importdata('Projections.mat');
angles      = importdata('Angles.mat');
model       = importdata('Local_classification_coord_OriOri.mat');
atoms       = importdata('Local_classification_type.mat');

% process the data: has to be double type in lsqcurvefit;
% the projections dimension has to be odd
projections = double(projections);
angles = angles';
projections = My_stripzero(projections,[280,280,51]);
projections_less = projections(2:end,2:end,:);

% set initial input parameters:
N = size(projections_less,1);   % calculate the dimension of projections
Res = 0.347;                    % indicate the pixel size for measured projections
CropHalfWidth = 4;              % the cropped bondary size for each atoms
% the atomic number for different type:
% use 28 for type 1, 45 for type 2, 78 for type 3
Z_arr = [28,45,78];

% running the least square fit for H and B factor estimation, then use them 
% to perform position refinement based on gradient descent

% set initial guess for H and B factor estimation
x0 = [50, 100, 150;...
    9, 9, 9];

% set the values into struct for lsqcurvefit
xdata = [];
xdata.model   = model;
xdata.atoms   = atoms;
xdata.angles  = angles;
xdata.Res     = Res;
xdata.Z_arr   = Z_arr;
xdata.volSize = N;
xdata.CropHalfWidth = CropHalfWidth;
xdata.halfWidth     = CropHalfWidth;
xdata.projections   = projections;

% set lower bound and upper bound
lb=[1, 1, 1;...
    1, 1, 1];
ub=[1000, 1000, 1000;...
    100, 100, 100];

opt = optimset('TolFun',1e-12);

% perform the least square fit for H and B factor estimation
[pXT, ~,~] = lsqcurvefit(@calcLinearPjs_fit, x0, xdata, projections_less, lb, ub, opt);
pXT = real(pXT);
disp(num2str(pXT))

% set the parameters for position refinement into structure
xdata.radius     = 1;       % radius of searching XYZ position, 1A is default
xdata.step_sz    = 20;      % step size can work up to 50, try it out
xdata.iterations = 10;      % iteration number
xdata.model_ori  = model;

% perform the position refinement
[y_pred,model_refined]    = gradient_XYZ_man(pXT, xdata, projections); 

model_refined_res = model_refined(:,:,end);
save('output/result.mat','model_refined_res');
