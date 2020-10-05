%% Main_2_match_mro_cell
% mro step 2£º get the final the mro network to fill in the whole nanoparticle 
% so that each mro solute atom corresponding at most 1 mro network, and make 
% the number of mro network as much as possible, the process begin from the 
% largest size mro network to the smallest

addpath('src/')
addpath('input/')

% read in files: finalized atomic types and the scaled boo parameter
atoms = importdata('final_atom_types.mat');
scaled_boo_para = importdata('scaled_boo_amorphous_region.mat');

% read in files: load the mro cell results calculate from mro step 1
load('mro_cell_amorphous_region_075_woTr.mat', 'mro_results','pos')
mro_results = mro_results(:,1:4);

% calculate the amorphous index corresponding to all atoms
cc = find(scaled_boo_para < 0.5 & atoms == 3);
cc_ind = zeros(numel(atoms),1);
cc_ind(scaled_boo_para < 0.5 & atoms == 3) = 1:numel(cc);

%% create lists: ind_array, connectivity_array (graph cell), rms_array
reverseStr = '';
% initialize the graph cell contained the connectivity for each mro cell
graph_cell = cell(size(mro_results,1),4);
for i = 1:size(graph_cell,1)
    for j = 1:4
        cluster_tmp = mro_results{i,j};
        cluster_graph_tmp = Cluster_graph(); % initialize the graph object
        % create the graph object for each mro cell
        if cluster_tmp(1).N_atoms>1
            cluster_tmp(1).grid_pts(:,1) = 0;
            temp_connectivity = cluster_tmp(1).neigh_matrix();
            cluster_graph_tmp.connectivity = cluster_tmp(1).neigh_matrix();
            cluster_graph_tmp.Graph = create_graph(cluster_graph_tmp.connectivity);
        end
        % create the graph object for all structure within each mro cell,
        % since for each type 3 atom, there could be more than 1 cell with
        % same largest network size
        if numel(cluster_tmp) > 1
            for k = 2:numel(cluster_tmp)
                cluster_graph_tmp_1 = Cluster_graph();
                if cluster_tmp(k).N_atoms>1
                cluster_tmp(k).grid_pts(:,1) = 0;
                cluster_graph_tmp_1.connectivity = cluster_tmp(k).neigh_matrix();
                cluster_graph_tmp_1.Graph = create_graph(cluster_graph_tmp_1.connectivity);
                end
                cluster_graph_tmp = [cluster_graph_tmp,cluster_graph_tmp_1];
            end
        end
        graph_cell{i,j} = cluster_graph_tmp;
    end
%     msg = sprintf(['Preparing graphic cell: atom num = ' sprintf('%i',i), '/' num2str(size(graph_cell,1))]);
    msg = sprintf('Preparing graphic cell: atom num = %i/%i', i, size(graph_cell,1));
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'),1,length(msg));
end
disp(' ')

% transfer the basic parameters from original mro cell to new mro graph cell; 
% calculate rms and shift for each mro graph cell
reverseStr = '';
for i = 1:size(graph_cell,1)
    for j = 1:size(graph_cell,2)
        cluster_graph_tmp = graph_cell{i,j};
        cluster_tmp = mro_results{i,j};
        for k = 1:numel(cluster_graph_tmp)
            cluster_tmp(k).grid_pts(:,1) = 0;
            cluster_graph_tmp(k).N_atoms = cluster_tmp(k).N_atoms;  % number of mro solute atoms
            cluster_graph_tmp(k).idxs = cluster_tmp(k).idxs;        % index for each solute atoms
            cluster_graph_tmp(k).type = cluster_tmp(k).type;        % type of mro network
            cluster_graph_tmp(k).param = cluster_tmp(k).param;      % lattice constant and rotation matrix
            cluster_graph_tmp(k).origin = cluster_tmp(k).origin;    % the center of mros
            cluster_graph_tmp(k).grid_pts = cluster_tmp(k).grid_pts;% the grid points for each solute atom in mro
            [shift,rms] = cluster_rms(cluster_tmp(k));              % calculate rms between each mro and template
            cluster_graph_tmp(k).shift = shift;
            cluster_graph_tmp(k).rms = rms;
        end
        graph_cell{i,j} = cluster_graph_tmp;
    end
    msg = sprintf('Calculating RMS: atom num = %i/%i', i, size(graph_cell,1));
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'),1,length(msg));
end
% for the solute atom which can produce more than 1 mro, only leave the one
% with smallest rms
for i = 1:size(graph_cell,1)
    for j = 1:size(graph_cell,2)
        cluster_graph_tmp = graph_cell{i,j};
        rms_arr = zeros(numel(cluster_graph_tmp),1);
        for k = 1:numel(cluster_graph_tmp)
            rms_arr(k) = cluster_graph_tmp(k).rms;
        end
        [~,rms_minInd] = min(rms_arr);
        graph_cell{i,j} = cluster_graph_tmp(rms_minInd);
    end
end
disp(' ')

N = 4; % single template match, 4 for all periodic template, 5 for all templates (include icosahedron)
% main while loop for splice the mro networks in the nanoparticle
% get the largest cluster
graph_cell  = graph_cell(:,1:4);
network_arr = {}; network_grid = {};
network_ind = []; network_type = [];
count = 1;
while true
    num_arr_arr = cellfun(@max_N_atoms,graph_cell);
    [dist_rms_arr,dist_rms_ind] = cellfun(@min_rms,graph_cell);
    maxVal = max(num_arr_arr(:));
    if maxVal <= 1
        break;
    end
    % get the index with largest mro size
    maxInd = find(num_arr_arr(:) == maxVal);
    % if there are more than 1 mro candidates, pick the one with smallest rms
    [~,minRMS_Ind] = min(dist_rms_arr(maxInd));
    [I1,I2] = ind2sub(size(num_arr_arr),maxInd(minRMS_Ind));
    I3 = dist_rms_ind(I1,I2);
    % reduce all the index in each group
    % reduce all connectivity and leave the largest cluster
    reduced_ind = graph_cell{I1,I2}(I3).idxs;
    if numel(reduced_ind) ~= maxVal
        keyboard;
    end
    network_arr  = [network_arr,reduced_ind];
    network_ind  = [network_ind;reduced_ind];
    network_type = [network_type,graph_cell{I1,I2}(I3).type];
    network_grid = [network_grid,graph_cell{I1,I2}(I3)];
    graph_cell  = cellfun(@(x)reduce_ind_group(x,reduced_ind),graph_cell,'UniformOutput',false);
    % record the type, template and vars
    % go back to while loop
    disp([num2str(count),':  Num: ',num2str(numel(reduced_ind)),'  Type: ',num2str(network_type(end))])
    count = count + 1;
end
num_arr_arr = cellfun(@numel,network_arr);

save('mro_list_amorphous_region_075.mat','network_arr','network_grid','network_ind','network_type','num_arr_arr');
