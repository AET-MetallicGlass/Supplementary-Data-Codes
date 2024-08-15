% HtBf_Par
clear
clc

pj_filename              = 'MG_multislice_refpj.mat';
angle_filename           = 'MG_multislice_refangle.mat';
atomtype_filename        = 'LocalC_MG_multislice_localtype.mat';
model_filename           = 'LocalC_MG_multislice_OriOri.mat';

pjs = 1:55;
model = importdata(model_filename);
atoms = importdata(atomtype_filename);

angles = importdata(angle_filename);
projections = importdata(pj_filename);
projections = max(projections,0);
projections = projections(2:end,2:end,pjs);

projections = My_paddzero(projections,size(projections)+[50 50 0],'double');

[N1,N2,num_pj] = size(projections);
halfWidth = 4;
Z_arr = [28 45 78];
Res = 0.347;


xdata = [];
xdata.Res     = Res;
xdata.Z_arr   = Z_arr;
xdata.halfWidth = halfWidth;
xdata.atoms = atoms;
xdata.model = model;
xdata.angles = angles;

%%
para0 = [1, 1.88, 3.96;...
    5.68, 5.52, 5.34];

model_refined = model;

lb=[1, 1, 1;...
    3, 3, 3];
ub=[1, 3, 5;...
    10, 10, 10];

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
    save([outputDir 'mg_ms_gs150_Fit_better_ini_HB_' num2str(jjjj) '.mat'],'y_pred','para0')
    
    xdata.projections=[];
    
    xdata.step_sz    = 1;
    xdata.iterations = 10;
    [y_pred,para0,errR] = gradient_B_2type_difB(para0, xdata, projections);
    save([outputDir 'mg_ms_gs150_Fit_better_ini_HB_itr10_' num2str(jjjj) '.mat'],'y_pred','para0','errR')
    
    xdata.step_sz    = 1;
    xdata.iterations = 10;
    [y_pred,para,errR] = gradient_fixHB_XYZ(para0, xdata, projections); model_refined = para(3:5,:);
    save([outputDir 'mg_ms_gs150_Fit_better_ini_HB_xyz_itr10_' num2str(jjjj) '.mat'],'y_pred','para','errR')
end
