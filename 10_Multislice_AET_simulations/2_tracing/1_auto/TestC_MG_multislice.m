
tracing_filename = 'polyn_tracing_MG_multislice_rad3_result%i_result.mat';
load('finalRec_MG_multislice_intp3.mat')
save_filename = 'peak_info_MG_multislice_rad%i_shift_result.mat';

classify_info = struct('Num_species', 3,  'halfSize',  3,  'plothalfSize',  1, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  0,  'separate_part',  120);

for i = 3
    temp_tracing_filename = sprintf(tracing_filename,i);
    load(temp_tracing_filename);
    atom_pos = TotPosArr(exitFlagArr==0,:)';
    
    b1 = find(atom_pos(1,:)<15 | atom_pos(1,:)>size(FinalVol_single,1)-15);
    b2 = find(atom_pos(2,:)<15 | atom_pos(2,:)>size(FinalVol_single,2)-15);
    b3 = find(atom_pos(3,:)<15 | atom_pos(3,:)>size(FinalVol_single,3)-15);
    
    bT = union(union(b1,b2),b3);
    atom_pos(:,bT) = [];
    
    [temp_model, temp_atomtype] = initial_class_kmean(...
        FinalVol_single, atom_pos, classify_info);
    
    [peak_info,intensity_arr,intensity_plot_arr] = ...
        plot_class_hist(FinalVol_single, temp_model, temp_atomtype, classify_info);
    
    temp_model_o = temp_model / 3 - 2 ;
    
    temp_save_filename = sprintf(save_filename,i);
    save(temp_save_filename,'temp_model_o',...
        'temp_atomtype','peak_info','intensity_arr','intensity_plot_arr')
end