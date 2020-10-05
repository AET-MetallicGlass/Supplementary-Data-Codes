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

% running the least square fit for H and B factor estimation, then use them 
% to perform position refinement based on gradient descent

% process the data: has to be double type in lsqcurvefit;
% the projections dimension has to be odd
projections = max(projections,0);
projections = projections(2:end,2:end,:);
projections = My_paddzero(projections,size(projections) + [50 50 0],'double');

[N1,N2,num_pj] = size(projections);
% the cropped bondary size for each atoms
halfWidth = 4;
% the atomic number for different type:
% use 28 for type 1, 45 for type 2, 78 for type 3
Z_arr   = [28 45 78];    
% indicate the pixel size for measured projections
Res     = 0.347;            

xdata = [];
xdata.Res       = Res;
xdata.Z_arr     = Z_arr;
xdata.halfWidth = halfWidth;
xdata.atoms     = atoms;
xdata.model     = model;
xdata.angles    = angles;

%%
% set initial guess for H and B factor estimation
para0 = [1, 1.36, 2.58;  13.3, 13.3, 13.3];
lb = [1, 1, 1;  5,  5,  5];
ub = [1, 2, 3;  15, 15, 15];

model_refined = model;

opt = optimset('TolFun',1e-12);

for jjjj=1:10
    fprintf('iteration num: %d; \n',jjjj);
    x0 = para0;
    x0(1,:)=x0(1,:)/x0(1,1);
    
    xdata.model = model_refined;
    xdata.model_ori = model_refined;
    xdata.projections=projections;
    
    [para0, ~,~] = lsqcurvefit(@Cal_Bproj_2type2, x0, xdata, projections, lb, ub, opt);
    fprintf('H1 = %.03f, H2 = %.03f, H3 = %.03f\n B11 = %.03f, B2 = %.03f, B3 = %.03f\n',...
        para0(1),para0(2),para0(3),para0(4),para0(5),para0(6));
    [y_pred,~] = Cal_Bproj_2type(para0, xdata, projections);
    
    xdata.projections = [];
    
    xdata.step_sz    = 1;
    xdata.iterations = 10;
    [y_pred,para0,errR] = gradient_B_2type_difB(para0, xdata, projections);
    
    xdata.step_sz    = 1;
    xdata.iterations = 10;
    [y_pred,para,errR] = gradient_fixHB_XYZ(para0, xdata, projections); model_refined = para(3:5,:);

end
model_refined_res = model_refined;
save('output/model_refined_res.mat','model_refined_res');
