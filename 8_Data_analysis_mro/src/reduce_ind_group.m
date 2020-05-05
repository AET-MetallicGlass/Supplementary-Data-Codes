function [cluster_graph_group] = reduce_ind_group(cluster_graph_group,reduced_ind)
for i = 1:numel(cluster_graph_group)
    cluster_graph_group(i) = reduce_ind(cluster_graph_group(i),reduced_ind);
end
end