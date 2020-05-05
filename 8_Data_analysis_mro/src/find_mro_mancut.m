function results = find_mro_mancut(idx,mro_type,pos,neigh,cut,TRUNCATE_SEARCH)
if nargin < 6
    TRUNCATE_SEARCH = 1;
end
q = [Cluster_mancut_new(pos(idx,:),idx,mro_type,cut)];
results = [q(1)];
combo_idx = 1;
visited   = cell(size(pos,1),1);
visited{idx} = [combo_idx];

while(size(q,2)>0)
    curr = q(1);
    q(1) = [];
    nn = [];
    for j = 1:size(curr.idxs,1)
        nn = [nn neigh{curr.idxs(j)}];
    end
    nn = unique(nn);
    
    if(TRUNCATE_SEARCH && rand >= min(1,100/size(q,1)))
        nn = [];
    end
    
    if(size(nn,1)*size(nn,2)==0)
        nn = [];
    end
    
    for j = 1:size(nn,2)
        if(~ismember(nn(j),curr.idxs))
            new_atom_idxs = sort([curr.idxs; nn(j)]);
            if(~checked_combo(visited, new_atom_idxs))
                combo_idx = combo_idx + 1;
                visited = add_combo(visited, new_atom_idxs, combo_idx);
                [child, code] = curr.add_pos(pos(nn(j),:),nn(j));
                if(code==1)
                    q = [q; child];
                    if(child.N_atoms == results(1).N_atoms)
                        results = [results; child];
                    elseif(child.N_atoms > results(1).N_atoms)
                        results = [child];
                    end
                end
            end
        end
    end
end
end

function checked = checked_combo(visited, idxs)
visited_sizes = zeros(size(idxs,1),1);
for i=1:size(visited_sizes,1)
    visited_sizes(i) = size(visited{idxs(i)},1);
end
[~,order] = sort(visited_sizes);
a = visited{idxs(order(1))};
if(size(a,1)*size(a,2) == 0)
    checked = 0;
    return
end
for i=2:size(idxs,1)
    a = intersect(a,visited{idxs(order(i))});
    if(size(a,1)*size(a,2) == 0)
        checked = 0;
        return
    end
end
checked = 1;
end

function visited = add_combo(visited, idxs, combo_idx)
for i=1:size(idxs,1)
    visited{idxs(i)} = [visited{idxs(i)}; combo_idx];
end
end

