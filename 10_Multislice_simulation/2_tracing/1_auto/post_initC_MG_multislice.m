load('peak_info_MG_multislice_rad3_shift_result.mat')
load('polyn_tracing_MG_multislice_rad3_result.mat')
%%
atom_pos = TotPosArr(exitFlagArr==0,:)';
atom_pos_all = atom_pos./3 - 2;

atom_pos_o = temp_model_o;
mean(atom_pos_o,2)
%%
load('Factor_MG_multislice.mat', 'final_Rec')
support_para = struct('th_dis_r_afterav', 0.90,   'dilate_size',  15, ...
                      'erode_size',       15,     'bw_size',      50000);

tight_support1 = obtain_tight_support(final_Rec,support_para);
support_para.erode_size = 25;
tight_support2 = obtain_tight_support(final_Rec,support_para);

%%
bdl_1 = 8;
bdl_2 = size(final_Rec,1) - 8;

ind_out1 = atom_pos_o(1,:) <= bdl_1 | atom_pos_o(2,:) <= bdl_1 | atom_pos_o(3,:) <= bdl_1;
ind_out2 = atom_pos_o(1,:) >= bdl_2 | atom_pos_o(2,:) >= bdl_2 | atom_pos_o(3,:) >= bdl_2;
atom_pos_o(:,ind_out1|ind_out2) = [];

ind_out1 = atom_pos_all(1,:) <= bdl_1 | atom_pos_all(2,:) <= bdl_1 | atom_pos_all(3,:) <= bdl_1;
ind_out2 = atom_pos_all(1,:) >= bdl_2 | atom_pos_all(2,:) >= bdl_2 | atom_pos_all(3,:) >= bdl_2;
atom_pos_all(:,ind_out1|ind_out2) = [];

%%
temp_pos_arr1 = [];
ind_arr1 = [];
for i = 1:size(atom_pos_o,2)
    temp_pos = round(atom_pos_o(:,i));
    if tight_support1(temp_pos(1),temp_pos(2),temp_pos(3)) == 1
        temp_pos_arr1 = [temp_pos_arr1, atom_pos_o(:,i)];
        ind_arr1 = [ind_arr1,i];
    end
end

ind_arr2 = [];
for i = 1:size(atom_pos_all,2)
    if min(pdist2(atom_pos_all(:,i)',temp_pos_arr1')) < 1e-4
        ind_arr2 = [ind_arr2,i];
    else
        temp_pos = round(atom_pos_all(:,i));
        if tight_support1(temp_pos(1),temp_pos(2),temp_pos(3)) == 1
            ind_arr2 = [ind_arr2,i];
        end
    end
end
%%
temp_pos_arr2 = atom_pos_all(:,ind_arr2);
save('post_initC_mg_yk_ms_gs150_ref2ar_175_pos_peakclass_arr_201003.mat','temp_pos_arr2')
