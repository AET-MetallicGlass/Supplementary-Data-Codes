%% croppedOut %%
% crops out m x m from the center of larger N x N array

function ROI = croppedOut(largeArray,cropSize)

n = size(largeArray);
nc = round((n+1)/2);

if length(cropSize) == 1
    cropSize = repmat(cropSize,length(n),1);
end
for ii = 1:length(n)
    vec = 1:cropSize(ii);
    cropC = round((cropSize(ii)+1)/2);
    cropVec{ii} = single(vec - cropC + nc(ii));
end

if length(n) == 2
    ROI = largeArray(cropVec{1}, cropVec{2});
elseif length(n) == 3
    ROI = largeArray(cropVec{1}, cropVec{2}, cropVec{3});
end

end
    