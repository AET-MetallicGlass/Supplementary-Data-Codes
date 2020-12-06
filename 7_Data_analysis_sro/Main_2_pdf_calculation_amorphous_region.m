%% Main_2_pdf_calculation_amorphous_region
% calculate both the pair distribution functions and partial pair distribution functions
% with amorphous atoms only inside the metallic glass nanoparticle by the scaled boo

addpath('src/')
addpath('input/')

% read in files: finalized atomic coordinates in Angstrom and types, as
% well as the scaled boo parameter
model = importdata('final_atom_coord_inA.mat');
atoms = importdata('final_atom_types.mat');
scaled_boo_amor = importdata('scaled_boo_amorphous_region.mat');
model = double(model);

% calculate the amorphous atoms and their types
cc = scaled_boo_amor < 0.5;
model = model(:,cc);
atoms = atoms(:,cc);

% set the parameters: step size and the range for pdf and ppdf
res = 1;  step = 0.05;  cutoff = 20;

% calculate the alpha shape of the nanoparticle
submodel = model - repmat(min(model,[],2),[1,size(model,2)]) + ones(size(model));
submodel = double(submodel);
shp      = alphaShape(submodel(1,:)',submodel(2,:)',submodel(3,:)',4);

% initialize a spherical index for later intersection calculation
[vMat,~] = spheretri(5000);
nsample  = size(vMat,1);

% initialize the radial distance and g(r) array for pdf
radius_arr = 0:step:cutoff;
pdf_arr    = zeros([size(radius_arr,2),1]);
volume_arr = zeros([size(radius_arr,2),1]);
ppdf_arr    = zeros([size(radius_arr,2),3,3]);

% perform the main part for pdf and ppdf calculation
for i=1:size(submodel,2)
    disp(num2str(i));
    
    temp_atom_arr   = atoms;
    temp_label      = temp_atom_arr(i);
    pdf_arr_temp    = zeros([size(radius_arr,2),1]);
    volume_arr_temp = zeros([size(radius_arr,2),1]);
    ppdf_arr_temp    = zeros([size(radius_arr,2),3,3]);
    
    for j = 1:size(submodel,2)
        
        dis = submodel(:,i) - submodel(:,j);
        dis = norm(dis,2);
        % for each pair, calculate the pair distribution and accumulate the
        % number to corresponding array postion
        if dis < cutoff
            ind = ceil(dis/step+0.01);
            pdf_arr_temp(ind) = pdf_arr_temp(ind)+1;
            ppdf_arr_temp(ind,temp_label,temp_atom_arr(j)) = ...
                ppdf_arr_temp(ind,temp_label,temp_atom_arr(j)) + 1;
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
    
    pdf_arr = pdf_arr + pdf_arr_temp;
    ppdf_arr = ppdf_arr + ppdf_arr_temp;
    volume_arr = volume_arr + volume_arr_temp;
    
end
% plot pdf as following code
shp = alphaShape(model(1,:)',model(2,:)',model(3,:)',ceil(4));
volume_alpha = volume(shp) + 4*pi()*90^2*0.1;
pdf_arr(1) = 0;
pdfnorm_arr= pdf_arr(:)./size(model,2)./(volume_arr(:)./volume_alpha);
radius_arr = radius_arr + 0.025;
pdf_amor05 = pdfnorm_arr;

pdf_amor05(1) = 0;
pdf_amor05(end) = pdf_amor05(end-1);
pdf_amor05 = imgaussfilt(pdf_amor05,2);
[base_g_sub1,ycorr_g_sub1] = baseline_normal(pdf_amor05);
pdf_amor05 = ycorr_g_sub1 + base_g_sub1/base_g_sub1(end);

figure(21); set(gcf,'position',[0,50,1000,600]); 
clf; hold on;
plot(r_dist_arr,pdf_amor05,'LineWidth',2);

