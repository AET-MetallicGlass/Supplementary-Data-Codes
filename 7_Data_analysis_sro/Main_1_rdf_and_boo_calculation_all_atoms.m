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
%% RDF plot
% plot rdf as following code
shp = alphaShape(model(1,:)',model(2,:)',model(3,:)',ceil(4));
volume_alpha = volume(shp) + 4*pi()*90^2*0.1;
rdf_arr(1) = 0;
rdfnorm_arr= rdf_arr(:)./size(model,2)./(volume_arr(:)./volume_alpha);
radius_arr = radius_arr + 0.05;
figure(11); plot(radius_arr,rdfnorm_arr);
%% BOO calculation
% the parameter for boo calculation is estimated by gaussian fitting of the
% 1st valley of the rdf of all atoms
Q_analysis_bondTh = 3.688;

% read in files: finalized atomic coordinates in pixel
model = importdata('final_atom_coord_inA.mat');

% perform Q4 and Q6 calculation, the periodic boundary condition is closed
% since it is nanoparticle
boo_Q4_full = obtain_Q_l_2ndNN_inds_PBC(model,Q_analysis_bondTh,4,[], 0);
boo_Q6_full = obtain_Q_l_2ndNN_inds_PBC(model,Q_analysis_bondTh,6,[], 0);

% the boo is calculated by scaling the square sum of Q4 and Q6 by the Q4
% and Q5 of perfect FCC lattice (calculating perfect FCC lattice give the Q4
% 0.190941 and Q6 0.574524)
scaled_SROP = sqrt(boo_Q4_full.^2 + boo_Q6_full.^2)./sqrt(0.190941^2 + 0.574524^2);

%% BOO plot
runsize = 35; falpha = 0.95;
edges{1} = linspace(0,0.4,runsize*2);
edges{2} = linspace(0,0.6,runsize);

fcc = [0.190941,0.574524,0.190941,0.574524];
bcc = [0.0363696,0.510688,0.0363696,0.510688];
hcp = [0.09722,0.484762,0.09722,0.0484762];

figure(12); clf; hold on; set(gcf,'position',[0,50,600,600]);

hh = histogram2(boo_Q4_full,boo_Q6_full,edges{1},edges{2},...
    'FaceColor','flat','EdgeColor','none');hold on

r1 = 0.5*norm(fcc(1:2),2);
r2 = 0.5*norm(fcc(1:2),2);
th = 0:pi/50:0.5*pi;
xunit = r1 * cos(th);
yunit = r2 * sin(th);
h1 = plot3(xunit, yunit,10000*ones(1,numel(yunit)),'--','Color','r','LineWidth',3);

r1 = 0.3*norm(fcc(1:2),2);
r2 = 0.3*norm(fcc(1:2),2);
th = 0:pi/50:0.5*pi;
xunit = r1 * cos(th);
yunit = r2 * sin(th);

xlim([-0.0 0.22]); ylim([0.0 0.60]); box on;

set(gca,'Xtick',0:0.1:0.2); set(gca,'XtickLabel',0:0.1:0.2);
set(gca,'Ytick',0:0.5:0.6); set(gca,'YtickLabel',0:0.5:0.6);
set(gca,'Ztick',0:50:100);  set(gca,'ZtickLabel',0:50:100)

axis square; view(2);
set(gca,'FontName','Arial','Fontsize',16,...
    'FontWeight','bold','linewidth',3,'TickLength',[0,0]);
shading interp;
ms = 100; grid off

scatter3(fcc(1),fcc(2),5,ms,'filled', 'MarkerEdgeColor','k','MarkerFaceColor','k');%[0.8500 0.3250 0.0980])
scatter3(bcc(1),bcc(2),5,ms,'filled','MarkerEdgeColor','k','MarkerFaceColor','k');%[0 0.4470 0.7410])
scatter3(hcp(1),hcp(2),5,ms,'filled','MarkerEdgeColor','k','MarkerFaceColor','k');%[0.4940 0.1840 0.5560])
material metal
shading interp
xlb = xlabel('Q4','FontSize',15,'FontName', 'Arial','FontAngle', 'normal','FontWeight','bold');
ylabel('Q6','FontSize',15,'FontName', 'Arial','FontWeight','bold')

text(fcc(1)-0.027,fcc(2),'fcc','FontSize',15,'FontName', 'Arial','FontWeight','bold')
text(bcc(1)+0.007,bcc(2),'bcc','FontSize',15,'FontName', 'Arial','FontWeight','bold')
text(hcp(1)+0.007,hcp(2),'hcp','FontSize',15,'FontName', 'Arial','FontWeight','bold')

h2 = colorbar;
ylabel(h2, 'Atom number','FontSize',15,'FontName', 'Arial','FontWeight','bold')
colormap('jet'); caxis([0,200]);

