%% Main_1_mro_cell_calculation
% mro step 1: calculate the potemtial mro based on breadth first search algorithm
% for each amorphous type 3 atom, calculate the largest potenitial type 3
% atoms network with specific lattice structure

addpath('src/')
addpath('input/')

% read in files: finalized atomic coordinates in Angstrom and types, as
% well as the scaled boo parameter
model = importdata('final_atom_coord_inA.mat');
atoms = importdata('final_atom_types.mat');
scaled_boo_para = importdata('scaled_boo_amorphous_region.mat');

% calculate the amorphous atoms and their types
cc = scaled_boo_para < 0.5 & atoms == 3;
model = model(:,cc);
atoms = atoms(:,cc);

% [Step 1] Load positions and calculate neighbor matrix  

R_UNIT  = 3.78;
r_inner = 3.78/R_UNIT;  % r_inner and r_outer are unitless 
r_outer = 6.09/R_UNIT; 
pos = model'/R_UNIT;

neigh = cell(size(pos,1),1);
for i = 1:size(pos,1)
    disp(num2str(i))
    for j = i + 1:size(pos,1)
        r_ij = pos(i,:) - pos(j,:);
        d_ij = sqrt(sum(r_ij.*r_ij));
        if(d_ij >= r_inner && d_ij <= r_outer)
            neigh{i} = [neigh{i},j];
            neigh{j} = [neigh{j},i];
        end
    end
end

% [Step 2] Calculate MRO for all atoms
cut = 0.75/R_UNIT;
mro_results = cell(size(pos,1),1);
for i = 1:size(pos,1)
    for j = 1:5
        mro_results{i} = find_mro_mancut(i,j,pos,neigh,cut,0);
    end
    disp(num2str(i))
end
