function [minVal,minInd] = min_rms(cluster_arr)
rms_arr = zeros(size(cluster_arr));
N_atoms_arr = zeros(size(cluster_arr));
for i = 1:numel(rms_arr)
    N_atoms_arr(i) = cluster_arr(i).N_atoms;
    rms_arr(i) = cluster_arr(i).rms;
end
max_N_atoms_id = find(N_atoms_arr == max(N_atoms_arr));
[minVal,minInd] = min(rms_arr(max_N_atoms_id));
minInd = max_N_atoms_id(minInd);
end