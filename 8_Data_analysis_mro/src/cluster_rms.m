function [clus_shift,clus_rms] = cluster_rms(obj)

% basis   = nn_basis(obj.type);
% a = obj.param(1);
% q = obj.param(2:5);
% rotm = quat_to_rotm(q);
% v = a*(basis')*(obj.grid_pts');
% r_grid = (rotm*v)';

[r_grid,~] = calc_r_grid(obj.param,obj.grid_pts,obj.type);
r_diff = obj.pos - obj.origin - r_grid;
clus_shift  = sqrt(sum(r_diff.*r_diff,2));

% [r_grid] = vars_to_grid_djc(cluster.param,cluster.type);
% r_grid_pos = cluster.grid_pts*r_grid + repmat(cluster.origin,[size(cluster.grid_pts,1),1]);
% clus_shift = sqrt(sum((r_grid_pos-cluster.pos).^2,2));

[clus_rms] = rms(clus_shift);
end