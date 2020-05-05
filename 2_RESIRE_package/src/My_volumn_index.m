%% My_volumn_index.m %%
% calculate corresponding index a small size matrix embedded in large size matrix

function vol_ind = My_volumn_index(big_size,ori_size)
vol_ind=zeros(length(big_size),2);

if length(big_size)~=length(ori_size)
    fprintf(1,'input volume dimension and and paddedsize length does not match!\n');
elseif sum(big_size < ori_size) > 0
    fprintf(1,'paddedsize should be equal to or smaller than original volume in all dimensions!\n');
else
    for i=1:length(big_size)
        if mod(big_size(i),2) == 0
            if mod(ori_size(i),2) == 0
                startind = 1 + (big_size(i)-ori_size(i))/2;
                endind = ori_size(i) + (big_size(i)-ori_size(i))/2;
            else
                startind = 1 + (big_size(i)-ori_size(i)+1)/2;
                endind = ori_size(i) + (big_size(i)-ori_size(i)+1)/2;
            end
        else
            if mod(ori_size(i),2) == 0
                startind = 1 + (big_size(i)-ori_size(i)-1)/2;
                endind = ori_size(i) + (big_size(i)-ori_size(i)-1)/2;
            else
                startind = 1 + (big_size(i)-ori_size(i))/2;
                endind = ori_size(i) + (big_size(i)-ori_size(i))/2;
            end
        end
        vol_ind(i,1) = startind;
        vol_ind(i,2) = endind;
    end
end

end