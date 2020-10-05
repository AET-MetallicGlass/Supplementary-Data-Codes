function [temp_model, temp_atomtype] = initial_class_kmean(rec, curr_model, classify_info)
% Yao Yang UCLA, 2020.08.20
% n-type classification by k-mean
% classify Non-atoms and atoms first
% then classify among real-atoms

if isfield(classify_info,'lnorm')            lnorm = classify_info.lnorm; %#ok<SEPEX>
else lnorm = 2;             %#ok<SEPEX>
end
if isfield(classify_info,'Num_species')      Num_species = classify_info.Num_species; %#ok<SEPEX>
else Num_species = 3;       %#ok<SEPEX>
end
if isfield(classify_info,'halfSize')         halfSize = classify_info.halfSize; %#ok<SEPEX>
else halfSize = 1;          %#ok<SEPEX>
end
if isfield(classify_info,'plothalfSize')     plothalfSize = classify_info.plothalfSize; %#ok<SEPEX>
else plothalfSize = 4;      %#ok<SEPEX>
end
if isfield(classify_info,'separate_part')    separate_part = classify_info.separate_part; %#ok<SEPEX>
else separate_part = 70;      %#ok<SEPEX>
end
if isfield(classify_info,'O_Ratio')          O_Ratio = classify_info.O_Ratio; %#ok<SEPEX>
else O_Ratio = 1;           %#ok<SEPEX>
end
if isfield(classify_info,'SPHyn')            SPHyn = classify_info.SPHyn; %#ok<SEPEX>
else SPHyn = 1;             %#ok<SEPEX>
end
if isfield(classify_info,'PLOT_YN')          PLOT_YN = classify_info.PLOT_YN; %#ok<SEPEX>
else PLOT_YN = 0;           %#ok<SEPEX>
end
%% generate points of intensities 

[box_inten] = get_box_intensity(rec, curr_model, halfSize, O_Ratio, SPHyn, 'linear');
[box_inten_plot] = get_box_intensity(rec, curr_model, plothalfSize, O_Ratio, SPHyn, 'linear');
%% k-mean main function
if lnorm == 2
    idx = kmeans(sum(box_inten,1)',Num_species+1,'Start','uniform');
%     idx = kmeans(points',Num_species,'Distance','sqeuclidean');
elseif lnorm == 1
    idx = kmeans(sum(box_inten,1)',Num_species+1,'Distance','cityblock','Start','uniform');
end
%% Alignment clustered type into correct species order
mean_arr = zeros(1,Num_species+1);
for i = 1:Num_species+1
    mean_arr(i) = mean(sum(box_inten(:,idx==i),1));
end
[~,sortMean] = sort(mean_arr);
for i = 1:Num_species+1
    idx(idx==sortMean(i)) = i+1000;
end
idx = idx-1000;
%% plot histogram
if PLOT_YN
    figure(203); clf;
    set(gcf,'Position',[0,0,400,900])
    % figure
    [hist_inten_plot,cen_integ_total_plot] = hist(sum(box_inten_plot,1),separate_part);
    y_up = round(max(hist_inten_plot)/10)*12;
    for i = 1:Num_species+2
        if i == 1
            subplot(Num_species+2,1,1);
            hist(sum(box_inten_plot,1),separate_part);
            title(sprintf('boxsize %d',halfSize*2+1));
        elseif i == 2
            intensity_integ_sub = sum(box_inten_plot(:,(idx==i-1)),1);
            subplot(Num_species+2,1,i)
            hist(intensity_integ_sub,cen_integ_total_plot);
            title(sprintf('%d Non-atoms',sum(idx==i-1)));
        else
            intensity_integ_sub = sum(box_inten_plot(:,(idx==i-1)),1);
            subplot(Num_species+2,1,i)
            hist(intensity_integ_sub,cen_integ_total_plot);
            title(sprintf('%d Type %d atoms',sum(idx==i-1),i-2));
        end
        xlabel('integrated intensity (a.u.)');
        ylabel('# atoms');
        
        ylim([0 y_up]);
        xlim([0 ceil(max(sum(box_inten_plot,1))/5)*5]);
    end
else
    fprintf('Rough classification:\n')
    fprintf('number of Non-atoms:  \t %d \n', sum(idx==1))
    for i = 1:Num_species
        fprintf('number of type%d atoms: \t %d \n', i, sum(idx==i+1))
    end
    fprintf('number of total atoms: \t %d \n', sum(idx~=1))
end
%% intermediate results
idx = idx - 1;

box_inten_type1 = get_box_intensity(rec, curr_model, halfSize, O_Ratio, 0, 'linear');
dim_b = round(size(box_inten_type1,1).^(1/3));
box_inten_type1 = reshape(box_inten_type1,[dim_b, dim_b, dim_b, numel(idx)]);
mean_box_type1  = mean(box_inten_type1(:,:,:,idx==1),4);

[atomtype, ~] = initial_class_L1norm(box_inten_type1, mean_box_type1, O_Ratio, halfSize, SPHyn );

temp_model          = curr_model(:,atomtype==1);
box_inten_sub       = box_inten(:,atomtype==1);
box_inten_plot_sub  = box_inten_plot(:,atomtype==1);

%% k-mean main function
if lnorm == 2
    idx = kmeans(sum(box_inten_sub,1)',Num_species,'Start','uniform');
%     idx = kmeans(points',Num_species,'Distance','sqeuclidean');
elseif lnorm == 1
    idx = kmeans(sum(box_inten_sub,1)',Num_species,'Distance','cityblock','Start','uniform');
end
%% Alignment clustered type into correct species order
mean_arr = zeros(1,Num_species);
for i = 1:Num_species
    mean_arr(i) = mean(sum(box_inten_sub(:,idx==i),1));
end
[~,sortMean] = sort(mean_arr);
for i = 1:Num_species
    idx(idx==sortMean(i)) = i+1000;
end
idx = idx-1000;
%% plot histogram
if PLOT_YN
    figure(204); clf;
    set(gcf,'Position',[0,0,400,900])
    % figure
    [hist_inten_plot_sub,cen_integ_total_plot_sub] = hist(sum(box_inten_plot_sub,1),separate_part);
    y_up = round(max(hist_inten_plot_sub)/10)*12;
    for i = 1:Num_species+1
        if i == 1
            subplot(Num_species+1,1,1);
            hist(sum(box_inten_plot_sub,1),separate_part);
            title(sprintf('boxsize %d',halfSize*2+1));
        else
            intensity_integ_sub = sum(box_inten_plot_sub(:,(idx==i-1)),1);
            subplot(Num_species+1,1,i)
            hist(intensity_integ_sub,cen_integ_total_plot_sub);
            title(sprintf('%d Type %d atoms',sum(idx==i-1),i-1));
        end
        xlabel('integrated intensity (a.u.)');
        ylabel('# atoms');
        
        ylim([0 y_up]);
        xlim([0 ceil(max(sum(box_inten_plot,1))/5)*5]);
    end
else
    fprintf('Initial classification:\n')
    for i = 1:Num_species
        fprintf('number of type%d atoms: \t %d \n', i, sum(idx==i))
    end   
    fprintf('number of total atoms: \t %d \n', numel(idx))
end
temp_atomtype = idx';
end