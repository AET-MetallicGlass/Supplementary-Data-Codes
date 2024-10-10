function [boxCoordinates, BoxCenter, BoxRadius, sphere] = create_box(BoxSize)
%M. Bartels, UCLA, 2014
%small function to create box coordinates and correpsonding spherical mask
%use odd BoxSize

BoxCenter = (BoxSize+1)/2;
BoxRadius = (BoxSize-1)/2;

%box coordinate systems
[boxX,boxY,boxZ]=ndgrid(-BoxRadius:BoxRadius,-BoxRadius:BoxRadius,-BoxRadius:BoxRadius);

%calculate spherical mask
sphere = sqrt(boxX.^2+boxY.^2+boxZ.^2)<=BoxSize/2;

boxCoordinates.x = boxX;
boxCoordinates.y = boxY;
boxCoordinates.z = boxZ;
end
