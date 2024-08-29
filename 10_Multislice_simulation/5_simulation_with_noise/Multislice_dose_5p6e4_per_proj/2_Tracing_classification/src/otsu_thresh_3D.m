function [ot,x]=otsu_thresh_3D(im)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Otsu Thresholding
% [ot,outim]=otsu_thresh(im)
%   where im is the grey image (integer 0-255)
%       Otsu Thresholding is perfromed on the image
%   ot is the threshold determined
%   outim is the image after thresholding

[dim1, dim2, dim3]=size(im);
pixels=dim1*dim2*dim3;

bins=[0:255]+0.5;
N=hist(reshape(im,pixels,1),bins);

Nnorm=N/sum(N); %normalising the bin frequencies to make probabilities
theta=cumsum(Nnorm); %cumulative probability
mu=cumsum(Nnorm.*[0:255]);

sigB2=(mu-mu(256)*theta).^2./(theta.*(1-theta)); %evaluate sigB2 over the t range

[p ot]=max(sigB2); %find the maximum value and the index where it is (this is the Otsu threshold)

x=(im>ot); %thresholding