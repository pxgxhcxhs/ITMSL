function [glacier_mask] = make_glacier_mask(dem,outline)
glacier_mask_z = zeros(size(dem.Z));
xx_min=min(outline.X);
xx_max=max(outline.X);
yy_min=min(outline.Y);
yy_max=max(outline.Y);

X=dem.georef.SpatialRef.XWorldLimits(1,1)+0.5*dem.georef.SpatialRef.CellExtentInWorldX(1,1):dem.georef.SpatialRef.CellExtentInWorldX(1,1):dem.georef.SpatialRef.XWorldLimits(1,2)-0.5*dem.georef.SpatialRef.CellExtentInWorldX(1,1);
Y=dem.georef.SpatialRef.YWorldLimits(1,1)+0.5*dem.georef.SpatialRef.CellExtentInWorldY(1,1):dem.georef.SpatialRef.CellExtentInWorldY(1,1):dem.georef.SpatialRef.YWorldLimits(1,2)-0.5*dem.georef.SpatialRef.CellExtentInWorldY(1,1);

[jj,ii]=size(glacier_mask_z);

for j=1:jj
    for i=1:ii
        if X(i)>=xx_min && X(i)<=xx_max && Y(j)>=yy_min && Y(j)<=yy_max
            in =inpolygon(X(i),Y(j),outline.X,outline.Y);
            if in==1
                glacier_mask_z(j,i) = 1;
            end
        end
    end
end
glacier_mask_z = flipud(glacier_mask_z);
glacier_mask = dem;
glacier_mask.Z = glacier_mask_z;
glacier_mask.name = 'glacier_mask';
end