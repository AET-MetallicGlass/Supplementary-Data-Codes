clear all
%%
Mag_shift=1;
shift=[0 0 0]';
GT_atoms=load('MG_final_atom_coord_inA.mat'); % Ground truth
Type=load('MG_final_atom_types.mat');
Type=Type.final_atom_types;
ID=find(Type==3);% Cal. RMSD of specified type of atoms. 1(Co), 2(Pd) 3(Pt)

N=length(ID);

ref_pos0 = GT_atoms.final_atom_coord_inAngstrom; 
ref_pos = ref_pos0(:,ID); % Cal. RMSD of different type of atoms.
%ref_pos = ref_pos0(:,:);% Cal. RMSD of all atoms.

ref_pos(1,:)=ref_pos(1,:)-mean(ref_pos(1,:));
ref_pos(2,:)=ref_pos(2,:)-mean(ref_pos(2,:));
ref_pos(3,:)=ref_pos(3,:)-mean(ref_pos(3,:));

%% Classified atoms
Model=load('Model_CoPdPt_Dose1p7e4.mat');
atoms=Model.temp_model;
atoms_type=Model.temp_atomtype;
ID2=find(atoms_type==3);% Cal. RMSD of specified type of atoms. 1(Co), 2(Pd) 3(Pt)

%%
atom_pos_all = atoms*0.347;% 
atom_pos=atom_pos_all(:,ID2); % Cal. RMSD of different type of atoms.
%atom_pos=atom_pos_all(:,1:end);% Cal. RMSD of all atoms.

atom_pos(1,:)=atom_pos(1,:)-mean(atom_pos(1,:));
atom_pos(2,:)=atom_pos(2,:)-mean(atom_pos(2,:));
atom_pos(3,:)=atom_pos(3,:)-mean(atom_pos(3,:));

%%
while Mag_shift>1e-4
m1 = ref_pos - shift;
m2= atom_pos(:,1:end);

difarr = [ ];
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
disp(num2str(numel(count_arr1)/size(m1,2)))

%%
figure
scatter3(ref_pos(1,:)-shift(1),ref_pos(2,:)-shift(2),ref_pos(3,:)-shift(3),50,'r.')
hold on
scatter3(atom_pos(1,:),atom_pos(2,:),atom_pos(3,:),20,'k.')

%%
dif=difarr;
sig=std(dif,0,2);

rmsd_value = sqrt(sum(sig.^2));
fprintf('rmsd = %.2fpm\n',rmsd_value*100);

figure(11); clf; set(gcf,'Position',[0,100,1000,450]);
hh = histogram(pdist2(double(dif'),[0,0,0]),0:0.005:0.8);
hh.FaceColor = [0.4 0.7 1];
hh.FaceAlpha = 1; 
xlim([0,0.80])

xlabel(['Deviation(',char(197),')'])
ylabel('Number of atoms')
set(gca,'FontSize',15,'FontName', 'Arial','FontWeight','bold','linewidth',2);
%% Output RMSD of different type of atoms by set ID= 1, 2, and 3 or output RMSD for all atoms
% save('dif_Dose1p7e4_AtomTypeAll.mat','dif')
% save('dif_Dose1p7e4_AtomType1.mat','dif')
% save('dif_Dose1p7e4_AtomType2.mat','dif')
% save('dif_Dose1p7e4_AtomType3.mat','dif')
return
%%
dif_AtomTypeAll=load('dif_Dose1p7e4_AtomTypeAll.mat');
dif_AtomTypeAll=dif_AtomTypeAll.dif;
dif_AtomType1=load('dif_Dose1p7e4_AtomType1.mat');
dif_AtomType1=dif_AtomType1.dif;
dif_AtomType2=load('dif_Dose1p7e4_AtomType2.mat');
dif_AtomType2=dif_AtomType2.dif;
dif_AtomType3=load('dif_Dose1p7e4_AtomType3.mat');
dif_AtomType3=dif_AtomType3.dif;


figure(11); clf; set(gcf,'Position',[0,0,350,900]);
subplot(3+1,1,1);
hh = histogram(pdist2(double(dif_AtomTypeAll'),[0,0,0]),0:0.01:0.8);
%hh = histogram(pdist2(double(dif(1,:)'),[0]),0:0.002:0.2);
hh.FaceColor = [0.4 0.7 1];
hh.FaceAlpha = 1;
xlim([0,0.80])
xlabel(['Deviation (',char(197),')'])
ylabel('Number of atoms')

subplot(3+1,1,2);
hh = histogram(pdist2(double(dif_AtomType1'),[0,0,0]),0:0.01:0.8);
%hh = histogram(pdist2(double(dif(1,:)'),[0]),0:0.002:0.2);
hh.FaceColor = [0.4 0.7 1];
hh.FaceAlpha = 1;
xlim([0,0.80])
xlabel(['Deviation (',char(197),')'])
ylabel('Number of atoms')

subplot(3+1,1,3);
hh = histogram(pdist2(double(dif_AtomType2'),[0,0,0]),0:0.01:0.8);

hh.FaceColor = [0.4 0.7 1];
hh.FaceAlpha = 1; 
xlim([0,0.80])
xlabel(['Deviation (',char(197),')'])
ylabel('Number of atoms')

subplot(3+1,1,4);
hh = histogram(pdist2(double(dif_AtomType2'),[0,0,0]),0:0.01:0.8);

hh.FaceColor = [0.4 0.7 1];
hh.FaceAlpha = 1; 
xlim([0,0.80])
xlabel(['Deviation (',char(197),')'])
ylabel('Number of atoms')

%title(sprintf('RMSD = 5 pm, Precision = 99.94%%'));


