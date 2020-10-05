%%
load('output/mro_list_amorphous_region_075.mat')

num_arr_larger5 = sum(num_arr_arr>=5);
hist_info = zeros(max(num_arr_arr),4);

for i = 1:max(num_arr_arr)
%     ind = num_arr_arr == i;
    for j = 1:4
    hist_info(i,j) = sum(num_arr_arr == i & network_type == j);
    end
end

figure
set(gcf,'Position',[0,0,1000,600])
hist_info(1:4,:) = 0;
c_Top = [...
    0.4940,    0.1840,    0.5560;...
    0.1,       0.4940,    0.0250;...
    0,         0.4470,    0.7410;...
    0.6350,    0.0780,    0.1840];
c_Top_1 = c_Top([3,4,2,1],:);

b = bar(hist_info(:,[3,4,2,1]),'stacked');

for k = 1:size(hist_info,2)
    b(k).FaceColor = c_Top_1(k,:);
    
end
xlim([4,24])
ylim([0,43])
xticks(5:2:25)
yticks(0:10:110)

xlabel('MRO-network size','FontSize',15,'FontName', 'Arial','FontWeight','bold')
ylabel('MRO-network number','FontSize',15,'FontName', 'Arial','FontWeight','bold')
box on
grid off
set(gca,'FontSize',15,'FontName', 'Arial','FontWeight','bold','linewidth',3);

%%
stat_tot_atom_info = zeros(1,4);
stat_tot_netw_info = zeros(1,4);
for i = 1:4
    stat_tot_atom_info(i) = sum(num_arr_arr(num_arr_arr>=5 & network_type == i));
    stat_tot_netw_info(i) = sum(num_arr_arr>=5 & network_type == i);
end

stat_tot_atom_info = stat_tot_atom_info([3,4,2,1]);
stat_tot_netw_info = stat_tot_netw_info([3,4,2,1]);
figure
set(gcf,'Position',[0,0,400,600])
clf
hold on

b = bar(stat_tot_atom_info./sum(stat_tot_atom_info),0.5,'FaceColor','flat');
for i = 1:4
b.CData(i,:) = c_Top_1(i,:);
end
xlim([0,5])
xticks(1:4)
xticklabels({'FCC','HCP','BCC','SC'})
ylim([0,0.42])
yticks(0:0.1:0.4)
% yticklabels({'0','1%','2%','3%'}) 
ylabel('MRO types ratio')
box on
set(gca,'FontSize',15,'FontName', 'Arial','FontWeight','bold','linewidth',3);
