% get_box_intensity
% Yao Yang UCLA 2020.08.20
function [points] = get_box_intensity(rec, curr_model, halfSize, O_Ratio, SPHyn, interp_type)

Num_atom = size(curr_model,2);

% obtain global intensity histogram
ds = 1/O_Ratio;

%return
[XX,YY,ZZ] = ndgrid(-halfSize: ds :halfSize, ...
                    -halfSize: ds :halfSize, ...
                    -halfSize: ds :halfSize);

%length(XX)

if SPHyn
    useInd = find( (XX.^2 + YY.^2 + ZZ.^2) <= (halfSize + 0.5*ds)^2 );
else
    useInd = 1:length(XX(:));
end

% size(XX)
% 
% length(XX(:))

% generate points coordinates
YY = YY(useInd); XX = XX(useInd); ZZ = ZZ(useInd);
y_set = zeros(length(XX(:)), Num_atom);
x_set = zeros(length(YY(:)), Num_atom);
z_set = zeros(length(ZZ(:)), Num_atom);

% size(y_set)
% interpolations for points
for k = 1:Num_atom
    y_set(:,k) = YY + curr_model(2,k);
    x_set(:,k) = XX + curr_model(1,k);
    z_set(:,k) = ZZ + curr_model(3,k);
end

% if strcmp(interp_type,'linear')
%     points = splinterp3(rec, y_set, x_set, z_set);
% else
points = interp3(rec, y_set, x_set, z_set, interp_type);
% end
%size(points)
end