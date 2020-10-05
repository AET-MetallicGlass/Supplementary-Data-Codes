function [temp_model, temp_atomtype] = local_class_kmean_sub(rec, curr_model, curr_types, classify_info)
% Yao Yang UCLA, 2020.08.05
% n-type local classification by k-mean
% only classify among real-atoms

if isfield(classify_info,'lnorm')            lnorm = classify_info.lnorm; %#ok<SEPEX>
else lnorm = 2;             %#ok<SEPEX>
end
if isfield(classify_info,'StopCri')          StopCri = classify_info.StopCri; %#ok<SEPEX>
else StopCri = 5;       %#ok<SEPEX>
end
if isfield(classify_info,'halfSize')         halfSize = classify_info.halfSize; %#ok<SEPEX>
else halfSize = 1;          %#ok<SEPEX>
end
if isfield(classify_info,'O_Ratio')          O_Ratio = classify_info.O_Ratio; %#ok<SEPEX>
else O_Ratio = 1;           %#ok<SEPEX>
end
if isfield(classify_info,'Radius')           Radius = classify_info.Radius; %#ok<SEPEX>
else Radius = 15;             %#ok<SEPEX>
end
if isfield(classify_info,'SPHyn')            SPHyn = classify_info.SPHyn; %#ok<SEPEX>
else SPHyn = 1;             %#ok<SEPEX>
end

%% generate points of intensities 

[box_inten] = get_box_intensity(rec, curr_model, halfSize, O_Ratio, SPHyn, 'linear');

num_types = numel(unique(curr_types));
%% main loop

endFlag = 0;  currDesc = [];
pre_atomtype = curr_types;
new_atomtype = zeros(size(curr_types));

while ~endFlag
fprintf('new round: \n');

for i = 1:size(curr_model,2)
    curr_atompos = curr_model(:,i);
    Dist = pdist2(curr_model',curr_atompos');
    BallInd = Dist ~= 0 & Dist < Radius;
    
    R_arr = zeros(num_types,1);
    for j = 1:num_types
        temp_type = pre_atomtype == j;
        R_temp_type = norm((box_inten(:,i) - mean(box_inten(:,BallInd' & temp_type),2)), lnorm);
%         R_temp_type = norm((sum(box_inten(:,i),1) - mean(box_inten(:,BallInd' & temp_type),1)), 1);
        R_arr(j) = R_temp_type;
    end
    
    [~, MinInd] = min(R_arr);
    new_atomtype(i) = MinInd;
end

for i = 1:num_types
    print_arr = ['num', num2str(i), ': %d; '];
    fprintf(print_arr,sum(new_atomtype == i));
end
fprintf('\n');

if sum(pre_atomtype ~= new_atomtype) == 0
    endFlag = 1;
    currDesc(end+1) = 0;
    pre_atomtype = new_atomtype;
else
    fprintf('discrepency: %d\n',sum(pre_atomtype ~= new_atomtype))
    currDesc(end+1) = sum(pre_atomtype~=new_atomtype);
    pre_atomtype = new_atomtype;
    if length(currDesc) > StopCri
        cutCri = currDesc(end-StopCri+1:end);
        if sum(cutCri == currDesc(end)) == length(cutCri)
            endFlag = 1;
        end
    end
end

end

temp_model = curr_model;
temp_atomtype = new_atomtype;
end