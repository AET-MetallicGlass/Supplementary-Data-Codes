
function [r_grid] = vars_to_grid_djc(vars,type)
a = vars(1);
q = vars(2:5);
rotm = quat_to_rotm(q);
v1 = [a,0,0];
v2 = [0,a,0];
v3 = [0,0,a];
u = rotm*[v1' v2' v3'];
basis = nn_basis(type);
r_grid = zeros(size(basis));
for i=1:size(r_grid,1)
    r_grid(i,:) = (u*(basis(i,:)'))';
end
end