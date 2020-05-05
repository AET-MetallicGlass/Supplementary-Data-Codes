%% Main_3_voronoi_calculation_amorphous_region
% calculate voronoi index for all atoms inside the nanoparticle

addpath('src/')
addpath('input/')

% read in files: finalized atomic coordinates in Angstrom and types, as
model = importdata('final_atom_coord_inA.mat');
model = model';

% set parameters: the cutoff area for the voronoi calculation to be 1% to
% reduce the thermal expansion effect
areacutoff = 0.01;
% the bondlength cutoff is achieved by gaussian fitting from the 1st valley
% of rdf curve of amorphous atoms in the nanoparticle
bondlengthcutoff = 3.78;

% perform the voronoi calcualtion
[V,R] = voronoin(model);

DT = delaunayTriangulation(model(:,1),model(:,2),model(:,3));
ConnectivityList = DT.ConnectivityList;

vor = []; vorid = []; vornk = []; vorarea = []; neigh = [];
for tid = 1:size(model,1)
    disp(num2str(tid))
    XR10 = V(R{tid},:);
    if sum(sum(isinf(XR10),1),2) >= 1
        ind=sum(ConnectivityList==tid,2)~=0;
        ind=ConnectivityList(ind,:);
        ind=unique(ind(:));
        dis=vecnorm(model(ind,:)-model(tid,:),2,2);
        ind=ind(dis<bondlengthcutoff);
        neigh{tid}=model(ind,:);
        continue
    end
    
    K = convhull(XR10);
    
    vec1 = XR10(K(:,1),:)-XR10(K(:,2),:);
    vec2 = XR10(K(:,1),:)-XR10(K(:,3),:);
    nk = cross(vec1,vec2);
    nk = nk./repmat(sqrt(sum(nk.^2,2)),[1 3]);
    
    ind = 1:size(nk,1);
    num = 0;
    faceid = [];
    facevor = [];
    while ~isempty(ind)
        flag = sum(nk(ind,:).*nk(ind(1),:),2);
        faceidtemp = abs(flag)>1-1e-5 & abs(sum((XR10(K(ind,1),:)-XR10(K(ind(1),1),:)).*nk(ind(1),:),2))<1e-5;
        tempid = K(ind(faceidtemp),:);
        num = num+1;
        faceid{num} = unique(tempid(:));
        
        % sort vertices of each face to follow clockwise or counterclockwise order
        % for plotting
        
        center = mean(XR10(faceid{num},:),1);
        pol = XR10(faceid{num},:)-center;
        pol = pol./repmat(sqrt(sum(pol.^2,2)),[1 3]);
        npol = size(pol,1);
        Y = dot(cross(repmat(pol(1,:),[npol,1]),pol)',repmat(nk(ind(1),:),[npol,1])');
        
        D = atan2d(Y,dot(repmat(pol(1,:),[npol,1])',pol'));
        D(D<0) = 360 + D(D<0);
        [~,sid] = sort(D);
        faceid{num} = faceid{num}(sid);
        ind(faceidtemp) = [];
    end
    facenk = []; facearea = [];
    for i=1:size(faceid,2)
        % calculate surface normal
        vec1 = XR10(faceid{i}(1),:)-XR10(faceid{i}(2),:);
        vec2 = XR10(faceid{i}(1),:)-XR10(faceid{i}(3),:);
        nk = cross(vec1,vec2);
        facenk(i,:) = nk/sqrt(sum(nk.^2));
        
        % calculate face area
        vec = XR10(faceid{i}(2:end),:)-XR10(faceid{i}(1),:);
        facearea(i) = 0.5*sum(vecnorm(cross(vec(1:end-1,:),vec(2:end,:)),2,2));
        facevor{i} = XR10(faceid{i},:);
    end
    vorid{tid}   = faceid;
    vor{tid}     = facevor;
    vornk{tid}   = facenk;
    vorarea{tid} = facearea;
    
    % remove face with small area
    faceremove = find(facearea < areacutoff*sum(facearea));
    ind = sum(ConnectivityList==tid,2)~=0;
    ind = ConnectivityList(ind,:);
    ind = unique(ind(:));
    
    for i = faceremove
        vec = XR10(faceid{i}(1),:)-model(tid,:);
        vec = sum(vec.*facenk(i,:))*facenk(i,:)*2;
        pos = [model(tid,:)+vec;model(tid,:)-vec];
        dis1 = vecnorm(model(ind,:)-pos(1,:),2,2);
        dis2 = vecnorm(model(ind,:)-pos(2,:),2,2);
        atomremove = find(dis1<1e-5 | dis2<1e-5);
        ind(atomremove) = [];
    end
    dis = vecnorm(model(ind,:)-model(tid,:),2,2);
    ind = ind(dis<bondlengthcutoff); % comment this line if don't want bond length regulation to the voronoi
    neigh{tid} = model(ind,:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % re-calculate voronoi after regulation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if length(ind)<5
        vor{tid}=[];
        continue
    end
    tidsub = find(ind==tid);
    [Vsub,Rsub] = voronoin(model(ind,:));
    XR10sub = Vsub(Rsub{tidsub},:);
    if sum(sum(isinf(XR10sub),1),2) >= 1
        vor{tid}=[];
        continue
    end
    Ksub = convhull(XR10sub);
    vec1 = XR10sub(Ksub(:,1),:)-XR10sub(Ksub(:,2),:);
    vec2 = XR10sub(Ksub(:,1),:)-XR10sub(Ksub(:,3),:);
    nk = cross(vec1,vec2);
    nk = nk./repmat(sqrt(sum(nk.^2,2)),[1 3]);
    
    indsub = 1:size(nk,1);
    num = 0;
    faceid = [];
    facevor = [];
    while ~isempty(indsub)
        flag = sum(nk(indsub,:).*nk(indsub(1),:),2);
        faceidtemp = abs(flag)>1-1e-5 & abs(sum((XR10sub(Ksub(indsub,1),:)-XR10sub(Ksub(indsub(1),1),:)).*nk(indsub(1),:),2))<1e-5;
        tempid = Ksub(indsub(faceidtemp),:);
        num = num+1;
        faceid{num} = unique(tempid(:));
        
        % sort vertices of each face to follow clockwise or counterclockwise order
        % for plotting
        
        center = mean(XR10sub(faceid{num},:),1);
        pol = XR10sub(faceid{num},:)-center;
        pol = pol./repmat(sqrt(sum(pol.^2,2)),[1 3]);
        npol = size(pol,1);
        Y = dot(cross(repmat(pol(1,:),[npol,1]),pol)',repmat(nk(indsub(1),:),[npol,1])');
        
        D = atan2d(Y,dot(repmat(pol(1,:),[npol,1])',pol'));
        D(D<0) = 360 + D(D<0);
        [~,sid] = sort(D);
        faceid{num} = faceid{num}(sid);
        indsub(faceidtemp) = [];
        facevor{num} = XR10sub(faceid{num},:);
    end
    vor{tid} = facevor;
    neigh{tid} = model(ind,:);
    
end

% calculate indices
% defined by voronoi volume

edgelengthcutoff = 0.00; % fraction of averaged edge length in a voronoi cell

indlist = zeros(1,6); % voronoi index list
Nindlist = 0; % counts of voronoi index
VoronoiID = [];
VoronoiID{1} = 0; % voronoi cell IDs in each index
badID = [];
facecounts = 0; % every voronoi cell has the same weight, not every face

voro_ind_atoms = zeros(size(model,1),6);
for tid = 1:size(model,1)
    facevor = vor{tid};
    if isempty(facevor)
        continue
    end
    edgelength = [];
    nedge = 0;
    totedgelen = 0;
    for i=1:size(facevor,2)
        ntemp = size(facevor{i},1);
        nedge = nedge + ntemp;
        edgelength{i} = vecnorm(facevor{i}-facevor{i}([ntemp 1:ntemp-1],:),2,2);
        if sum(edgelength{i}>6)
            badID = [badID tid];
            continue
        end
        totedgelen = totedgelen + sum(edgelength{i});
    end
    ind = zeros(1,6);
    for i = 1:size(facevor,2)
        n = sum(edgelength{i} >= totedgelen / nedge * edgelengthcutoff);
        if n <= 6
            ind(n) = ind(n)+1;
        end
    end
    facecounts = facecounts + ind/size(neigh{tid},1);
    
    voro_ind_atoms(tid,:) = ind;
    temp = indlist == ind;
    id = find(sum(temp,2) == 6);
    if ~isempty(id)
        Nindlist(id) = Nindlist(id)+1;
        VoronoiID{id} = [VoronoiID{id} tid];
    else
        indlist = [indlist;ind];
        Nindlist = [Nindlist,1];
        VoronoiID{end+1} = [tid];
    end
end