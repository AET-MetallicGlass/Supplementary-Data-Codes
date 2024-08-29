inputDir = '';
pj_filename              = [inputDir 'MG_multislice_refpj.mat'];
angle_filename           = [inputDir 'MG_multislice_refangle.mat'];
support_filename         = [inputDir 'MG_multislice_supp.mat'];

outputDir = '';
results_filename         = [outputDir 'GD_MG_multislice.mat'];
filenameFinalRecon       = [outputDir 'Factor_MG_multislice.mat'];

%% GraDIRE parameters

GraDIRE = GraDIRE_Reconstructor();

%%% See the README for description of parameters

GraDIRE.filename_Projections = pj_filename;
GraDIRE.filename_Angles = angle_filename ;
GraDIRE.filename_Support = support_filename;
GraDIRE.filename_Results = results_filename;

GraDIRE.oversamplingRatio = 4;
GraDIRE.numIterations = 200; 
GraDIRE.monitor_R = 1;
GraDIRE.monitorR_loopLength = 20;

GraDIRE.griddingMethod = 1; 

%GENFIRE.vector1 = [0 0 1];
%GENFIRE.vector2 = [0 1 0];
GraDIRE.vector3 = [1 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin GENFIRE

GraDIRE = readFiles(GraDIRE);
% GraDIRE.InputProjections = My_paddzero(GraDIRE.InputProjections,[320,320,size(GraDIRE.InputProjections,3)]);
GraDIRE = CheckPrepareData(GraDIRE);
GraDIRE = runGridding(GraDIRE); 
GraDIRE = reconstruct(GraDIRE);
% SaveResults(GENFIRE);
%%


[Rfactor,Rarray,simuProjs] = calc_Rfactor_realspace_general(...
    GraDIRE.reconstruction,...
    GraDIRE.InputProjections,...
    GraDIRE.InputAngles,...
    [0 0 1],[0 1 0],[1 0 0]);
final_Rec = GraDIRE.reconstruction;
final_Rfactor = Rfactor;
final_Rarray = Rarray;
final_simuProjs = simuProjs;

save(filenameFinalRecon, 'final_Rec','final_Rfactor','final_Rarray','final_simuProjs');
% save(results_filename, 'GraDIRE','-v7.3');

