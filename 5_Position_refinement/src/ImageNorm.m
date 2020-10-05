function FinalImage = ImageNorm(Image,Value)
FinalImage = zeros(size(Image));
if numel(size(Image)) == 3
    for i = 1:size(Image,3)
        sumImage=sum(sum(Image(:,:,i)));
        FinalImage(:,:,i) = Image(:,:,i) * Value / sumImage;
    end
else
    if numel(size(Image)) == 2
        sumImage=sum(sum(Image));
        FinalImage = Image * Value / sumImage;
    else
        fprintf('Image can only be 2D or 3D! /%');
    end
end
end