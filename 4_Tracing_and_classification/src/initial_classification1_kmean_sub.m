function [temp_model, temp_atomtype] = initial_classification1_kmean_sub(...
    RecVol_padded, curr_model, classify_info)
% n-type classification by k-mean
% only classify among real-atoms

if isfield(classify_info,'lnorm')            lnorm = classify_info.lnorm; %#ok<SEPEX>
else lnorm = 2;             %#ok<SEPEX>
end
if isfield(classify_info,'Num_species')      Num_species = classify_info.Num_species; %#ok<SEPEX>
else Num_species = 3;       %#ok<SEPEX>
end
if isfield(classify_info,'phase_shift_flag') phase_shift_flag = classify_info.phase_shift_flag; %#ok<SEPEX>
else phase_shift_flag = 1;  %#ok<SEPEX>
end
if isfield(classify_info,'halfSize')         halfSize = classify_info.halfSize; %#ok<SEPEX>
else halfSize = 1;          %#ok<SEPEX>
end
if isfield(classify_info,'plothalfSize')     plothalfSize = classify_info.plothalfSize; %#ok<SEPEX>
else plothalfSize = 4;      %#ok<SEPEX>
end
if isfield(classify_info,'SPHyn')            SPHyn = classify_info.SPHyn; %#ok<SEPEX>
else SPHyn = 1;             %#ok<SEPEX>
end
if isfield(classify_info,'separate_part')    separate_part = classify_info.separate_part; %#ok<SEPEX>
else separate_part = 100;   %#ok<SEPEX>
end
if isfield(classify_info,'PLOT_YN')          PLOT_YN = classify_info.PLOT_YN; %#ok<SEPEX>
else PLOT_YN = 0;           %#ok<SEPEX>
end
%% obtain global intensity histogram

[xX,yY,zZ] = ndgrid(-halfSize:halfSize,-halfSize:halfSize,-halfSize:halfSize);
SphereInd = find(xX.^2+yY.^2+zZ.^2 <=(halfSize+0.5)^2);
[xXp,yYp,zZp] = ndgrid(-plothalfSize:plothalfSize,-plothalfSize:plothalfSize,-plothalfSize:plothalfSize);
SphereInd_plot = find(xXp.^2+yYp.^2+zZp.^2 <=(plothalfSize+0.5)^2);

if SPHyn
  useInd = SphereInd;
  useInd_plot =SphereInd_plot;
else
  useInd = 1:length(xX);
  useInd_plot = 1:length(xXp);
end

%% compute initial average Fe atom, Pt atom and non atom

% re-produce array for integrated intensity
intensity_integ_plot = zeros(1,size(curr_model,2));

% integrate intensity for each traced peak
for j=1:size(curr_model,2)
    curr_pos  = round(curr_model(:,j));
    box_integ = RecVol_padded(curr_pos(1)-plothalfSize:curr_pos(1)+plothalfSize,...
        curr_pos(2)-plothalfSize:curr_pos(2)+plothalfSize,curr_pos(3)-plothalfSize:curr_pos(3)+plothalfSize);
    intensity_integ_plot(j) = sum(box_integ(useInd_plot));
end

%% generate points of intensities 
Num_atom = size(curr_model,2);
pt_size  = length(useInd);
points   = zeros(pt_size, Num_atom);
for k = 1:size(curr_model,2)
    curr_x = round(curr_model(1,k));
    curr_y = round(curr_model(2,k));
    curr_z = round(curr_model(3,k));
    
    if phase_shift_flag == 1
        dx = curr_model(1,k) - curr_x;
        dy = curr_model(2,k) - curr_y;
        dz = curr_model(3,k) - curr_z;
    end
    curr_vol = RecVol_padded(curr_x-halfSize-1:curr_x+halfSize+1,...
        curr_y-halfSize-1:curr_y+halfSize+1,...
        curr_z-halfSize-1:curr_z+halfSize+1);
    if phase_shift_flag == 1
        shiftModel = FourierShift3D_Yao(curr_vol,dx,dy,dz);
        shift_vol = shiftModel(2:end-1,2:end-1,2:end-1);
    else
        shift_vol = curr_vol(2:end-1,2:end-1,2:end-1);
    end
    
    points(:,k) = shift_vol(useInd);
end
%% k-mean main function
if lnorm == 2
    idx = kmeans(points',Num_species);
%     idx = kmeans(points',Num_species,'Distance','sqeuclidean');
elseif lnorm == 1
    idx = kmeans(points',Num_species,'Distance','cityblock');
end
%% Alignment clustered type into correct species order
mean_arr = zeros(1,Num_species);
for i = 1:Num_species
    mean_arr(i) = mean(sum(points(:,idx==i),1));
end
[~,sortMean] = sort(mean_arr);
for i = 1:Num_species
    idx(idx==sortMean(i)) = i+1000;
end
idx = idx-1000;
%% plot histogram
if PLOT_YN
figure(203); clf;
set(gcf,'Position',[0,0,400,900])
% figure
[hist_inten_plot,cen_integ_total_plot]= hist(intensity_integ_plot,separate_part);
y_up = round(max(hist_inten_plot)/10)*12;
for i = 1:Num_species+1
    if i == 1
        subplot(Num_species+1,1,1);
        hist(intensity_integ_plot,separate_part);
        title(sprintf('boxsize %d',halfSize*2+1));
    else
        
        intensity_integ_sub = intensity_integ_plot(idx==i-1);
        subplot(Num_species+1,1,i)
        hist(intensity_integ_sub,cen_integ_total_plot);
        title(sprintf('%d Type %i atoms',sum(idx==i-1),i-1));
    end
    xlabel('integrated intensity (a.u.)');
    ylabel('# atoms');
    
    ylim([0 y_up]);
    xlim([0 ceil(max(intensity_integ_plot)/5)*5]);
end
end
%% final results
temp_model = curr_model;
temp_atomtype = idx';
end