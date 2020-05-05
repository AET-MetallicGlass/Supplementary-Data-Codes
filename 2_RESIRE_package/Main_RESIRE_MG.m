%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                %%
%%                     Welcome to RESIRE!                         %%
%%          REal Space Iterative Reconstruction Engine            %%
%%                                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('src/')
addpath('src/splinterp/')

% add the path to load measured projections and angles, you can comment it 
% and move the projections and angles into input folder
addpath('../1_Measured_data/') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                          User Parameters                              %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% See the object RESIRE_Reconstructor() for description of parameters
% please run this on super computer due to the large array size and
% oversampling ratio

pj_filename         = '../1_Measured_data/Projections.mat';
angle_filename      = '../1_Measured_data/Angles.mat';

results_filename    = 'output/RESIRE_experiment_result.mat';

RESIRE = RESIRE_Reconstructor();

RESIRE.filename_Projections    = pj_filename;
RESIRE.filename_Angles         = angle_filename ;
RESIRE.filename_Results        = results_filename;

RESIRE = RESIRE.set_parameters(...
    'oversamplingRatio' ,4    ,'numIterations'       ,200 ,... 
    'monitor_R'         ,true ,'monitorR_loopLength' ,20 ,... 
    'griddingMethod'    ,1    ,'vector3'             ,[1 0 0], ...
    'use_parallel'      ,1);

RESIRE = readFiles(RESIRE);
RESIRE = CheckPrepareData(RESIRE);
RESIRE = runGridding(RESIRE); 
RESIRE = reconstruct(RESIRE);

RESIRE = ClearCalcVariables(RESIRE);

Reconstruction = RESIRE.reconstruction;

save('output/Reconstruction.mat','Reconstruction')
% SaveResults(RESIRE);
