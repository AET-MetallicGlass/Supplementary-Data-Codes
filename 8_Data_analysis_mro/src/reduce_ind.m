function [cluster_graph_tmp] = reduce_ind(cluster_graph_tmp,reduced_ind)
idxs = cluster_graph_tmp.idxs;
conn_arr = cluster_graph_tmp.connectivity;
G = cluster_graph_tmp.Graph;
grid_pts = cluster_graph_tmp.grid_pts;

if isempty(G)
    return;
end
if numel(idxs) == 1 || numel(conn_arr) == 1
    cluster_graph_tmp.N_atoms = 0;
    cluster_graph_tmp.idxs = [];
    cluster_graph_tmp.type = 0;
    return
end
[~,ia,~] = intersect(idxs,reduced_ind);
for i = 1:numel(ia)
    %     keyboard;
    % assume connectivity matrix conn_arr
    ind = find(conn_arr(i,:) == 1);
    G = rmedge(G,i,ind);
end

% connectivity = conn_arr;
% conn_arr(ia,:) = 0;
% conn_arr(:,ia) = 0;

% G = create_graph(conn_arr);

[bin] = conncomp(G);
bintype = sort(unique(bin));
binsize = zeros(size(bintype));

[~,sub_ind] = setdiff(idxs,reduced_ind);
for i = 1:numel(bintype)
    binsize(i) = sum(bin(sub_ind)==bintype(i));
end

[~,maxInd] = max(binsize);
graph_cluster_ind = sub_ind(bin(sub_ind) == bin(maxInd));
SG = subgraph(G, graph_cluster_ind);
SG = simplify(SG);
if ~isempty(intersect(graph_cluster_ind,ia))
    keyboard;
end
cluster_graph_tmp.N_atoms = numel(graph_cluster_ind);
cluster_graph_tmp.idxs = idxs(graph_cluster_ind);
cluster_graph_tmp.connectivity = conn_arr(graph_cluster_ind,graph_cluster_ind);
cluster_graph_tmp.Graph = SG;
cluster_graph_tmp.grid_pts = grid_pts(graph_cluster_ind,:);
if numel(graph_cluster_ind) == 1
    cluster_graph_tmp.shift = 0;
    cluster_graph_tmp.rms = 0;
else
    cluster_graph_tmp.shift = cluster_graph_tmp.shift(graph_cluster_ind);
    cluster_graph_tmp.rms = rms(cluster_graph_tmp.shift);
end
end