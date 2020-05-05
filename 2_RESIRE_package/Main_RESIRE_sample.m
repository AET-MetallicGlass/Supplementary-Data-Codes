%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                %%
%%                     Welcome to RESIRE!                         %%
%%          REal Space Iterative Reconstruction Engine            %%
%%                                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath('src/')
addpath('src/splinterp/')
% package splinterp is a fast C++ libary for parallel calculation of
% linear, bilinear and trilinear interpolation which need C++ compiler,
% please check https://github.com/apryor6/splinterp to find more detail

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                          User Parameters                              %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% See the object RESIRE_Reconstructor() for description of parameters

pj_filename         = 'input/sample_projections_amorphous.mat';
angle_filename      = 'input/sample_angles_amorphous.mat';

results_filename    = 'output/sample_amorphous_res.mat';

RESIRE = RESIRE_Reconstructor();

RESIRE.filename_Projections    = pj_filename;
RESIRE.filename_Angles         = angle_filename ;
RESIRE.filename_Results        = results_filename;

RESIRE = RESIRE.set_parameters(...
    'oversamplingRatio' ,3    ,'numIterations'       ,100 ,... 
    'monitor_R'         ,true ,'monitorR_loopLength' ,10 ,... 
    'griddingMethod'    ,1    ,'vector3'             ,[1 0 0], ...
    'use_parallel'      ,1);

RESIRE = readFiles(RESIRE);
RESIRE = CheckPrepareData(RESIRE);
RESIRE = runGridding(RESIRE); 
RESIRE = reconstruct(RESIRE);

RESIRE = ClearCalcVariables(RESIRE);

Reconstruction = RESIRE.reconstruction;

SaveResults(RESIRE);
