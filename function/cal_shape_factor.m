function [shape_factor,xxyy] = cal_shape_factor(glacier_mask,outline,initial_ice_thickness)
%Calculate glacier shape factor

%Relevant citations to read are:

%Ramsankaran, Raaj, A. Pandit, and M. F. Azam. "Spatially Distributed Ice-Thickness 
%Modelling for Chhota Shigri Glacier in Western Himalayas, India." Article, 
%International Journal of Remote Sensing 39, no. 10 (2018): 3320-43. 
%https://dx.doi.org/10.1080/01431161.2018.1441563.

%%%%%%%%%%%% Extracting glacial skeletons
BW = bwmorph(glacier_mask.Z,'skel',Inf);

%figure
%imshow(BW)
%%%%%%%%%%%% Calculate the coordinates of the center point of the skeleton pixel
[m,n]=find(BW==1);
for i=1:size(m)
    xx(i,1) = glacier_mask.georef.SpatialRef.XWorldLimits(1)+n(i)*glacier_mask.georef.SpatialRef.CellExtentInWorldX-1/2*glacier_mask.georef.SpatialRef.CellExtentInWorldX;
    yy(i,1) = glacier_mask.georef.SpatialRef.YWorldLimits(2)-m(i)*glacier_mask.georef.SpatialRef.CellExtentInWorldY+1/2*glacier_mask.georef.SpatialRef.CellExtentInWorldY;
    hh(i,1) = initial_ice_thickness(m(i),n(i));
    half_width(i,1) = dis2poly([xx(i,1); yy(i,1)],[(outline.X)' (outline.Y)']);
end

%%%%%%%%%%% Remove points outside the outline
[in, ~]=inpolygon(xx,yy,outline.X,outline.Y);
ind =find(in==0);
xx(ind,:)=[]; yy(ind,:)=[]; m(ind,:)=[]; n(ind,:)=[];

shape_factor = nanmean(2/pi*atan(half_width./hh));

xxyy = [xx yy];
end

