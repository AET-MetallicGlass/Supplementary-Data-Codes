function [avatom] = compute_average_atom_from_vol(DataMatrix,atom_pos,atom_ind,boxhalfsize)

% dbstop if isempty(atom_ind);
if isempty(atom_ind)
    avatom = zeros(boxhalfsize*2+1,boxhalfsize*2+1,boxhalfsize*2+1);
else
size_atom_ind = size(atom_ind);
if size_atom_ind(1) ~= 1
    atom_ind = atom_ind';
end
    
total_vol = zeros(boxhalfsize*2+1,boxhalfsize*2+1,boxhalfsize*2+1);

for kkk=atom_ind
    curr_x = round(atom_pos(1,kkk));
    curr_y = round(atom_pos(2,kkk));
    curr_z = round(atom_pos(3,kkk));
    
    curr_vol = DataMatrix(curr_x-boxhalfsize:curr_x+boxhalfsize,...
                            curr_y-boxhalfsize:curr_y+boxhalfsize,...
                            curr_z-boxhalfsize:curr_z+boxhalfsize);
    
    total_vol = total_vol + curr_vol;
end

    avatom = total_vol / length(atom_ind);
end
end