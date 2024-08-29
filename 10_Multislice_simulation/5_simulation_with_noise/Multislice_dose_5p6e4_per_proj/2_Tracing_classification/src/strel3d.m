function se=strel3d(sesize)

% function se=STREL3D(sesize)
%
% STREL3D creates a 3D sphere as a structuring element. Three-dimensional 
% structuring elements are much better for morphological reconstruction and
% operations of 3D datasets. Otherwise the traditional MATLAB "strel"
% function will only operate on a slice-by-slice approach. This function
% uses the aribtrary neighborhood for "strel."
% 
% Usage:        se=STREL3D(sesize)
%
% Arguments:    sesize - desired diameter size of a sphere (any positive 
%               integer)
%
% Returns:      the structuring element as a strel class (can be used
%               directly for imopen, imclose, imerode, etc)
% 
% Examples:     se=strel3d(1)
%               se=strel3d(2)
%               se=strel3d(5)
%
% 2014/09/26 - LX 
% 2014/09/27 - simplification by Jan Simon


sw=(sesize-1)/2; 
ses2=ceil(sesize/2);            % ceil sesize to handle odd diameters
[y,x,z]=meshgrid(-sw:sw,-sw:sw,-sw:sw); 
m=sqrt(x.^2 + y.^2 + z.^2); 
b=(m <= m(ses2,ses2,sesize)); 
se=strel('arbitrary',b);