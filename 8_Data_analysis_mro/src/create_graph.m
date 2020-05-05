function G = create_graph(connectivity)
G = graph;
G = addnode(G,size(connectivity,1));
if numel(connectivity) ~= 1
% for i = 1:size(connectivity,1)
%     G = addedge(G,i,find(connectivity(i,:) == 1));
% end
% G = simplify(G);

for i = 1:size(connectivity,1)
    ind = find(connectivity(i,:) == 1);
    for j = 1:numel(ind)
        try
            G = addedge(G,i,ind(j));
        catch
        end
    end
end
try
    G = simplify(G);
catch
end
end