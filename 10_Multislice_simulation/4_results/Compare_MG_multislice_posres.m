clear
Mag_shift=1;
shift=[0 0 0]';
load('mg_ms_gs150_Fit_better_ini_HB_xyz_itr10_10.mat', 'para')
ref_pos = para(3:5,:);
while Mag_shift>1e-4
m1 = ref_pos - shift;

m2=importdata('final_atom_coord_inA.mat');

difarr = [];
count_arr1 = [];
count_arr2 = [];
for i=1:size(m2,2)
    dif=m1-m2(:,i);
    dis=sqrt(sum(dif.^2,1));
    [dis,ind]=min(dis);
    if dis<=1.5
        difarr=[difarr dif(:,ind)];
        count_arr1 = [count_arr1,ind];
        count_arr2 = [count_arr2,i];
    end
end
shift=shift+mean(difarr,2);
Mag_shift=sum(abs(mean(difarr,2)));
end
disp(num2str(numel(count_arr1)/size(m2,2)))
%%
dif=difarr;
sig=std(dif,0,2);
rmsd_value = sqrt(sum(sig.^2));
fprintf('rmsd = %.2fpm\n',rmsd_value*100);
%%
figure(11); clf; set(gcf,'Position',[0,100,1000,450]);
hh = histogram(pdist2(double(dif'),[0,0,0]),0:0.02:1);
xlim([0,1]); xticks(0:0.2:1)
ylim([0,2250]); yticks(0:500:2000)
xlabel(['Deviation(',char(197),')'])
ylabel('Number of atoms')
%%
local_atomtype = importdata('final_atom_types.mat');

local_type_ms = importdata('LocalC_MG_multislice_localtype.mat');
local_type_ms_match = cell(3,1);
local_type_ms_match{3} = local_type_ms(count_arr1(local_atomtype(count_arr2)==3));
local_type_ms_match{2} = local_type_ms(count_arr1(local_atomtype(count_arr2)==2));
local_type_ms_match{1} = local_type_ms(count_arr1(local_atomtype(count_arr2)==1));
%%
load('initialC_shift_MG_multislice_rad3.mat')

figure(12); clf; set(gcf,'Position',[0,100,1000,450]); hold on;
b1 = bar(peak_info(1,:),peak_info(3,:));
b2 = bar(peak_info(1,:),peak_info(4,:));
b3 = bar(peak_info(1,:),peak_info(5,:));
xlabel('intensity (a.u.)');
ylabel('Number of atoms'); box on;
legend({'Type-1','Type-2','Type-3'})
xlim([0,4*10^4]);
%%
datainfo = zeros(5,3);
datainfo(1,1) = sum(cellfun(@(x)numel(x),local_type_ms_match));
datainfo(1,2) = 18356;
datainfo(2,1) = sum(local_type_ms_match{1}==1);
datainfo(2,2) = sum(local_atomtype==1);
datainfo(3,1) = sum(local_type_ms_match{2}==2);
datainfo(3,2) = sum(local_atomtype==2);
datainfo(4,1) = sum(local_type_ms_match{3}==3);
datainfo(4,2) = sum(local_atomtype==3);
datainfo(5,1) = sum(datainfo(2:4,1));
datainfo(5,2) = 18356;
datainfo(:,3) = datainfo(:,1)./datainfo(:,2);

fprintf('1: position matched atoms ratio = \t%.02f%%\n',datainfo(1,3)*100);
fprintf('2: matched type 1 atoms ratio = \t%.02f%%\n',datainfo(2,3)*100);
fprintf('3: matched type 2 atoms ratio = \t%.02f%%\n',datainfo(3,3)*100);
fprintf('4: matched type 3 atoms ratio = \t%.02f%%\n',datainfo(4,3)*100);
fprintf('5: matched all type atoms ratio = \t%.02f%%\n',datainfo(5,3)*100);
%%
figure(13); clf; set(gcf,'Position',[0,100,1000,450]);

for i = 1:3
    subplot(1,3,i);
    sig=std(dif(:,local_atomtype(count_arr2)==i),0,2);
    rmsd_value = sqrt(sum(sig.^2));
    fprintf('rmsd = %.2fpm\n',rmsd_value*100);

    hh = histogram(pdist2(double(dif(:,local_atomtype(count_arr2)==i)'),[0,0,0]),0:0.02:1);
    xlim([0,1]); xticks(0:0.2:1);
    ylim([0,1250]); yticks(0:500:2000);
    xlabel(['Deviation(',char(197),')'])
    ylabel('Number of atoms')
end
