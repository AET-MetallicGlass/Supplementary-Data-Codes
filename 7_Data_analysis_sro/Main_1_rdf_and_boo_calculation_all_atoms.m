%% Main_1_rdf_and_boo_calculation_all_atoms
% calculate the radial distribution functions and bond orientation order
% with all atoms inside the metallic glass nanoparticle

addpath('src/')
addpath('input/')

% read in files: finalized atomic coordinates in Angstrom and types
model = importdata('final_atom_coord_inA.mat');
atoms = importdata('final_atom_types.mat');

% set the parameters: step size and the range for rdf and pdf
step = 0.1;  cutoff = 20;

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

% the parameter for boo calculation is estimated by gaussian fitting of the
% 1st valley of the rdf of all atoms
Q_analysis_bondTh = 3.74;

% read in files: finalized atomic coordinates in pixel
model_ori = importdata('final_atom_model_inPixel.mat');
% the pixel size is 0.347A and the center of reconstruction volume is 161
model_ori = (model_ori-161)*0.347;

% perform Q4 and Q6 calculation, the periodic boundary condition is closed
% since it is nanoparticle
boo_Q4_full = obtain_Q_l_2ndNN_inds_PBC(model_ori,Q_analysis_bondTh,4,[], 0);
boo_Q6_full = obtain_Q_l_2ndNN_inds_PBC(model_ori,Q_analysis_bondTh,6,[], 0);

% the boo is calculated by scaling the square sum of Q4 and Q6 by the Q4
% and Q5 of perfect FCC lattice (calculating perfect FCC lattice give the Q4
% 0.190941 and Q6 0.574524)
scaled_SROP = sqrt(boo_Q4_full.^2 + boo_Q6_full.^2)./sqrt(0.190941^2 + 0.574524^2);
