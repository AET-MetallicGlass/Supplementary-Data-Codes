%% Main_2_rdf_calculation_amorphous_region
% calculate both the radial distribution functions and pair distribution functions
% with amorphous atoms only inside the metallic glass nanoparticle by the scaled boo

addpath('src/')
addpath('input/')

% read in files: finalized atomic coordinates in Angstrom and types, as
% well as the scaled boo parameter
model = importdata('final_atom_coord_inA.mat');
atoms = importdata('final_atom_types.mat');
scaled_boo_amor = importdata('scaled_boo_amorphous_region.mat');

% calculate the amorphous atoms and their types
cc = scaled_boo_amor < 0.5;
model = model(:,cc);
atoms = atoms(:,cc);

% set the parameters: step size and the range for rdf and pdf
res = 1;  step = 0.1;  cutoff = 30;

% calculate the alpha shape of the nanoparticle
submodel = model - repmat(min(model,[],2),[1,size(model,2)]) + ones(size(model));
submodel = double(submodel);
shp      = alphaShape(submodel(1,:)',submodel(2,:)',submodel(3,:)',4);

% initialize a spherical index for later intersection calculation
[vMat,~] = spheretri(5000);
nsample  = size(vMat,1);

% initialize the radial distance and g(r) array for rdf
radius_arr = 0:step:cutoff;
rdf_arr    = zeros([size(radius_arr,2),1]);
volume_arr = zeros([size(radius_arr,2),1]);
pdf_arr    = zeros([size(radius_arr,2),3,3]);

% perform the main part for rdf and pdf calculation
for i=1:size(submodel,2)
    disp(num2str(i));
    
    temp_atom_arr   = atoms;
    temp_label      = temp_atom_arr(i);
    rdf_arr_temp    = zeros([size(radius_arr,2),1]);
    volume_arr_temp = zeros([size(radius_arr,2),1]);
    pdf_arr_temp    = zeros([size(radius_arr,2),3,3]);
    
    for j = 1:size(submodel,2)
        
        dis = submodel(:,i) - submodel(:,j);
        dis = norm(dis,2);
        % for each pair, calculate the pair distribution and accumulate the
        % number to corresponding array postion
        if dis < cutoff
            ind = ceil(dis/step+0.01);
            rdf_arr_temp(ind) = rdf_arr_temp(ind)+1;
            pdf_arr_temp(ind,temp_label,temp_atom_arr(j)) = ...
                pdf_arr_temp(ind,temp_label,temp_atom_arr(j)) + 1;
        end
        
    end
    
    % calculate the intersection between the spherical volume with
    % nanoparticle as the scaled factor
    for j = 0:step:cutoff
        spoints = vMat*(j+step/2)+repmat(submodel(:,i)',[size(vMat,1),1]);
        in = inShape(shp,spoints(:,1),spoints(:,2),spoints(:,3));
        ind = round(j/step)+1;
        volume_arr_temp(ind) = sum(in) / nsample * ( 4/3 * pi() * ((j+step)^3-j^3) );
    end
    
    rdf_arr = rdf_arr + rdf_arr_temp;
    pdf_arr = pdf_arr + pdf_arr_temp;
    volume_arr = volume_arr + volume_arr_temp;
    
end
