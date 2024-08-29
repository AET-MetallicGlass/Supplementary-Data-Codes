% Polynomial tracing based on reconstruction of original orientation
clear all; close all
addpath('./src');
%% load data
RESIRE=load('Resire_CoPdPt_Dose1p7e4_rec.mat')
FinalVol=RESIRE.OBJ.reconstruction;

load("polyn_tracing_CoPdPt_Dose1p7e4_result.mat") % 
atom_pos0 = TotPosArr(exitFlagArr==0,:); %
%% initial speration of type-2 (Pd) and type-3 (Pt)

atom_pos1=atom_pos0(1:end,:)/3-2; 

classify_info1 = struct('Num_species', 3,  'halfSize',  1,  'plothalfSize',  3, ...
      'O_Ratio', 3, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  200, 'lnorm',2);

[temp_model1, temp_atomtype1] = initial_class_kmean_sub(FinalVol, atom_pos1', classify_info1);  
%% speration of type-1 (Co) and non-atoms
ID_type2=find(temp_atomtype1==2);
ID_type3=find(temp_atomtype1==3);
ID_type3_type2=[ID_type2 ID_type3];

atom_pos2=atom_pos1(:,:);
atom_pos2(ID_type3_type2,:)=[];
atom_pos_type3_type2=atom_pos1(ID_type3_type2,:);

classify_info2 = struct('Num_species', 2,  'halfSize',  1,  'plothalfSize',  3, ...
      'O_Ratio', 3, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  200, 'lnorm',2);

[temp_model2, temp_atomtype2] = initial_class_kmean_sub(FinalVol, atom_pos2', classify_info2);

ID_type1=find(temp_atomtype2==2);
atom_pos_type1=atom_pos2(ID_type1,:);
%% speration of type-1 (Co), type-2 (Pd) type-3 (Pt)
Atom_pos_all=[atom_pos_type3_type2; atom_pos_type1];

classify_info = struct('Num_species', 3,  'halfSize',  3,  'plothalfSize',  3, ...
      'O_Ratio', 3, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  200, 'lnorm',2);

[temp_model, temp_atomtype] = initial_class_kmean_sub(FinalVol, Atom_pos_all', classify_info);
return
%% Output
save('Model_CoPdPt_Dose1p7e4.mat',"temp_atomtype","temp_model");










