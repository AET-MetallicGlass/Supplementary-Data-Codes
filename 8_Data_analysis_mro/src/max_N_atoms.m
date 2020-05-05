function [maxVal,maxInd] = max_N_atoms(cluster_arr)
rms_arr = zeros(size(cluster_arr));
N_atoms_arr = zeros(size(cluster_arr));
for i = 1:numel(rms_arr)
    N_atoms_arr(i) = cluster_arr(i).N_atoms;
end
[maxVal,maxInd] = max(N_atoms_arr);
end