function [effect_k] = cal_effect_k(dem_basal,shape_factor,dem,glacier_mask)
%Calculate the effective stress coefficient

%% Search for the locations with the maximum and minimum elevation values.
[m1,n1] = find(dem.Z == nanmax(nanmax(dem.Z)));
[m2,n2] = find(dem.Z == nanmin(nanmin(dem.Z)));

dh = nanmax(nanmax(dem.Z))-nanmin(nanmin(dem.Z));
dist = (sqrt((m1-m2)^2+(n1-n2)^2))*dem.cellsize;
alpha = atand(dh/dist);

x1 = dem.georef.SpatialRef.XWorldLimits(1,1)+(n1+0.5)*dem.georef.SpatialRef.CellExtentInWorldX(1,1);
y1 = dem.georef.SpatialRef.YWorldLimits(1,2)-(m1+0.5)*dem.georef.SpatialRef.CellExtentInWorldY(1,1);

x2 = dem.georef.SpatialRef.XWorldLimits(1,1)+(n2+0.5)*dem.georef.SpatialRef.CellExtentInWorldX(1,1);
y2 = dem.georef.SpatialRef.YWorldLimits(1,2)-(m2+0.5)*dem.georef.SpatialRef.CellExtentInWorldY(1,1);

aspect_glacier = atand((y2-y1)/(x2-x1));
%% Determine the quadrant and calculate the azimuth.
if x2>x1 && y2>y1
    
elseif x2>x1 && y2<y1
    aspect_glacier = abs(aspect_glacier)+90;
elseif x2<x1 && y2<y1
    aspect_glacier = 270-aspect_glacier;
elseif x2<x1 && y2>y1
    aspect_glacier = 270-aspect_glacier;
end

%% Determine the aspect and extract the maximum upward angle of the flow direction.
aspect_basal = aspect(dem_basal);
basal_slope = gradient8(dem_basal,'deg');
if aspect_glacier>=0 && aspect_glacier<22.5  %%% N
    [m,n] = find(aspect_basal.Z>=157.5 & aspect_basal.Z<202.5);
    for i=1:length(m)
        aspect_basal.Z(m(i,1),n(i,1)) = 999;
    end
    aspect_basal.Z(aspect_basal.Z<999) = nan;
    msk_file = aspect_basal.Z./aspect_basal.Z;
    basal_slope.Z = basal_slope.Z.*msk_file;
    slope_max = nanmax(nanmax(basal_slope.Z));
elseif aspect_glacier>=337.5 && aspect_glacier<360 %%% N
    [m,n] = find(aspect_basal.Z>=157.5 & aspect_basal.Z<202.5);
    for i=1:length(m)
        aspect_basal.Z(m(i,1),n(i,1)) = 999;
    end
    aspect_basal.Z(aspect_basal.Z<999) = nan;
    msk_file = aspect_basal.Z./aspect_basal.Z;
    basal_slope.Z = basal_slope.Z.*msk_file;
    slope_max = nanmax(nanmax(basal_slope.Z));
elseif aspect_glacier>=22.5 && aspect_glacier<67.5 %%% NE
    [m,n] = find(aspect_basal.Z>=202.5 & aspect_basal.Z<247.5);
    for i=1:length(m)
        aspect_basal.Z(m(i,1),n(i,1)) = 999;
    end
    aspect_basal.Z(aspect_basal.Z<999) = nan;
    msk_file = aspect_basal.Z./aspect_basal.Z;
    basal_slope.Z = basal_slope.Z.*msk_file;
    slope_max = nanmax(nanmax(basal_slope.Z));
    
elseif aspect_glacier>=67.5 && aspect_glacier<112.5 %%% E
    [m,n] = find(aspect_basal.Z>=247.5 & aspect_basal.Z<292.5);
    for i=1:length(m)
        aspect_basal.Z(m(i,1),n(i,1)) = 999;
    end
    aspect_basal.Z(aspect_basal.Z<999) = nan;
    msk_file = aspect_basal.Z./aspect_basal.Z;
    basal_slope.Z = basal_slope.Z.*msk_file;
    slope_max = nanmax(nanmax(basal_slope.Z));
    
elseif aspect_glacier>=112.5 && aspect_glacier<157.5 %%% SE
    [m,n] = find(aspect_basal.Z>=292.5 & aspect_basal.Z<337.5);
    for i=1:length(m)
        aspect_basal.Z(m(i,1),n(i,1)) = 999;
    end
    aspect_basal.Z(aspect_basal.Z<999) = nan;
    msk_file = aspect_basal.Z./aspect_basal.Z;
    basal_slope.Z = basal_slope.Z.*msk_file;
    slope_max = nanmax(nanmax(basal_slope.Z));
elseif aspect_glacier>=157.5 && aspect_glacier<202.5 %%% S
    [m1,n1] = find(aspect_basal.Z>=337.5 & aspect_basal.Z<360);
    [m2,n2] = find(aspect_basal.Z>=0 & aspect_basal.Z<22.5);
    m = [m1;m2]; n = [n1;n2];
    for i=1:length(m)
        aspect_basal.Z(m(i,1),n(i,1)) = 999;
    end
    aspect_basal.Z(aspect_basal.Z<999) = nan;
    msk_file = aspect_basal.Z./aspect_basal.Z;
    basal_slope.Z = basal_slope.Z.*msk_file;
    slope_max = nanmax(nanmax(basal_slope.Z));    
elseif aspect_glacier>=202.5 && aspect_glacier<247.5 %%% SW
    [m,n] = find(aspect_basal.Z>=22.5 & aspect_basal.Z<67.5);
    for i=1:length(m)
        aspect_basal.Z(m(i,1),n(i,1)) = 999;
    end
    aspect_basal.Z(aspect_basal.Z<999) = nan;
    msk_file = aspect_basal.Z./aspect_basal.Z;
    basal_slope.Z = basal_slope.Z.*msk_file;
    slope_max = nanmax(nanmax(basal_slope.Z));
elseif aspect_glacier>=247.5 && aspect_glacier<292.5 %%% W
    [m,n] = find(aspect_basal.Z>=67.5 & aspect_basal.Z<112.5);
    for i=1:length(m)
        aspect_basal.Z(m(i,1),n(i,1)) = 999;
    end
    aspect_basal.Z(aspect_basal.Z<999) = nan;
    msk_file = aspect_basal.Z./aspect_basal.Z;
    basal_slope.Z = basal_slope.Z.*msk_file;
    slope_max = nanmax(nanmax(basal_slope.Z));
elseif aspect_glacier>=292.5 && aspect_glacier<337.5 %%% NW
    [m,n] = find(aspect_basal.Z>=112.5 & aspect_basal.Z<157.5);
    for i=1:length(m)
        aspect_basal.Z(m(i,1),n(i,1)) = 999;
    end
    aspect_basal.Z(aspect_basal.Z<999) = nan;
    msk_file = aspect_basal.Z./aspect_basal.Z;
    basal_slope.Z = basal_slope.Z.*msk_file;
    slope_max = nanmax(nanmax(basal_slope.Z));
end
beta = slope_max + alpha;

%% calculate effect_k  f*sin(surface_slope)/tan(beta)<=effect_k<1
surface_slope = gradient8(dem,'deg');
effect_k = shape_factor*sind(surface_slope.Z)/tand(beta);
effect_k = (effect_k+glacier_mask.Z)/2;
k = dem;
k.Z = effect_k;
effect_k = moving_average(k);
effect_k = effect_k.Z;
end

