
Recon_intp_filename     = 'finalRec_MG_multislice_intp3.mat';
filename_initialClass   = 'initialC_shift_MG_multislice_rad3.mat';
filename_localClass     = 'localC_shift_MG_multislice_rad3.mat';

newmodel        = importdata('finalmodel_MG_multislice.mat');
FinalVol_single = importdata(Recon_intp_filename);

classify_info = struct('Num_species', 3,  'halfSize',  3,  'plothalfSize',  1, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  0,  'separate_part',  70);

new_model_L = (newmodel -20 +2).*3;

[temp_model, temp_atomtype] = initial_class_kmean_sub(...
    FinalVol_single, new_model_L, classify_info);

[peak_info,intensity_arr,intensity_plot_arr] = ...
    plot_class_hist(FinalVol_single, temp_model, temp_atomtype, classify_info);

temp_model_o = temp_model / 3 - 2 ;

save(filename_initialClass,'temp_model_o',...
    'temp_atomtype','peak_info','intensity_arr','intensity_plot_arr')

classify_info.Radius = 10/0.347*3;
[temp_model, local_atomtype] = local_class_kmean_sub(FinalVol_single, temp_model, temp_atomtype, classify_info);

save(filename_localClass,'local_atomtype');
