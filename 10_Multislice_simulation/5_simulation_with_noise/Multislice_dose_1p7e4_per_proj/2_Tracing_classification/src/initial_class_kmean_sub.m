function [temp_model, temp_atomtype] = initial_class_kmean_sub(rec, curr_model, classify_info)
% Yao Yang UCLA, 2020.08.05
% n-type classification by k-mean
% only classify among real-atoms

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

box_inten_plot=box_inten_plot/216;%(2*halfSize)^3;
%% k-mean main function
if lnorm == 2
    idx = kmeans(sum(box_inten,1)',Num_species,'Start','uniform');
%     idx = kmeans(points',Num_species,'Distance','sqeuclidean');
elseif lnorm == 1
    idx = kmeans(sum(box_inten,1)',Num_species,'Distance','cityblock','Start','uniform');
end
%% Alignment clustered type into correct species order
mean_arr = zeros(1,Num_species);
for i = 1:Num_species
    mean_arr(i) = mean(sum(box_inten(:,idx==i),1));
end
[~,sortMean] = sort(mean_arr);
for i = 1:Num_species
    idx(idx==sortMean(i)) = i+1000;
end
idx = idx-1000;
%% plot histogram
if PLOT_YN
figure(203); clf;
%set(gcf,'Position',[0,0,400,900])
set(gcf,'Position',[0,0,350,900])
% figure
[hist_inten_plot,cen_integ_total_plot] = hist(sum(box_inten_plot,1),separate_part);
y_up = round(max(hist_inten_plot)/10)*12;
for i = 1:Num_species+1
    if i == 1
        subplot(Num_species+1,1,1);
        hist(sum(box_inten_plot,1),separate_part);
        title(sprintf('boxsize %d',halfSize*2+1));
h = findobj(gca,'Type','patch');
% h.FaceColor = [0 0.5 0.5];
% h.EdgeColor = 'w';

h.FaceColor = [0.4 0.7 1];
h.FaceAlpha = 1;

    else
        intensity_integ_sub = sum(box_inten_plot(:,(idx==i-1)),1);
        subplot(Num_species+1,1,i)
        hist(intensity_integ_sub,cen_integ_total_plot);
        title(sprintf('%d Type %i atoms',sum(idx==i-1),i-1));

h = findobj(gca,'Type','patch');
% h.FaceColor = [0 0.5 0.5];
% h.EdgeColor = 'w';

h.FaceColor = [0.4 0.7 1];
h.FaceAlpha = 1;
%set(gca,'FontSize',12,'FontName', 'Arial','FontWeight','bold','linewidth',2);
% xlim([0,0.80])
% %xticks(0:0.2:1)
% %ylim([0,100])
% %yticks(0:500:2000)
% xlabel(['Deviation(',char(197),')'])
% ylabel('Number of atoms')
% set(gca,'FontSize',15,'FontName', 'Arial','FontWeight','bold','linewidth',2);

    end
    xlabel('Integrated intensity (a.u.)');
  %  ylabel('# atoms');
    ylabel('Number of atoms');
    
    ylim([0 y_up]);
    xlim([0 ceil(max(sum(box_inten_plot,1))/5)*5]);
end
else
    fprintf('Initial classification sub:\n')
    for i = 1:Num_species
        fprintf('Number of type%d atoms: \t %d \n', i, sum(idx==i))
    end   
    fprintf('Number of total atoms: \t %d \n', numel(idx))
end
%% final results
temp_model = curr_model;
temp_atomtype = idx';
end