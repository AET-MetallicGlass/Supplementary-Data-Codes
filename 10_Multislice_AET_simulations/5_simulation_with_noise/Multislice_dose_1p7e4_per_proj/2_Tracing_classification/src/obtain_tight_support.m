%% Obtain tight support by parameter structure
% modified from Yongsoo's code My_obtain_tight_support_ver1()
% manually setting the otsu threshold, dilate size and erode_size
% can also remove small isolated area by bw_size
% Author: Yao Yang,  Latest update: 2020.08.05
% Coherent Imaging Group, UCLA

function curr_Supportt = obtain_tight_support(RECvol,para_info)

if isfield(para_info,'th_dis_r_afterav')    
    th_dis_r_afterav = para_info.th_dis_r_afterav;
else th_dis_r_afterav = 0.90;
end
if isfield(para_info,'dilate_size')    
    dilate_size = para_info.dilate_size;
else dilate_size = 11;
end
if isfield(para_info,'erode_size')    
    erode_size = para_info.erode_size;
else erode_size = 11;
end
fprintf('ostu_threshold: %0.2f, dilate size: %d, erode size: %d\n',th_dis_r_afterav,dilate_size,erode_size)
% smooth the volume
curr_RECvol = smooth3(RECvol,'b',9);

% Otsu threshold
im=double(curr_RECvol/max(curr_RECvol(:))*255);
[ot,~]=otsu_thresh_3D(im);
ot=ot*max(curr_RECvol(:))/255*th_dis_r_afterav;

curr_Support = (curr_RECvol>ot) * 1;

% make the mask slightly larger
se =strel3d(3);
curr_Support = imdilate(curr_Support,se);  

% make the mask quite larger
se =strel3d(dilate_size);
curr_Support = imdilate(curr_Support,se);  

if isfield(para_info,'bw_size')    
    bw_size = para_info.bw_size;
    curr_Support = bwareaopen(curr_Support,bw_size);
end

% make the mask quite smaller
se =strel3d(erode_size);
curr_Supportt = imerode(curr_Support,se);  
end