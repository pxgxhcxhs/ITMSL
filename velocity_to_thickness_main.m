%%%% Ice Thickness Model considering Sliding Law, main
clear
clc
%% Step 0 - set path
path_to_file = 'data\';
result_path = 'result\';

%% Step 1: Initialize parameters
n = 3;  %Glen's flow constant
Ac = 3.24e-24;  %Arrhenius creep constant Pa^(-3)*s^(-1)
f = 0.8; %lateral drag correction ('shape factor')
p_i = 900; %density of ice  kg/m^3
g = 9.8; %gravity m/s^2
slope_minimum_threhsold = 4;
slope_uncertainties = 5;
B = 465; %ice cerrp parameter 1/(MPa^3*yr)
cc = 0.852;

%% Step 2: Load the ice, topography, and velocity, temperature
%Load full dem
deeeeem = GRIDobj('data/dem.tif');

%load glacier outline file
outline=shaperead('data/outline.shp');

%make glacier mask file
glacier_mask = make_glacier_mask(deeeeem,outline);
dem_glacier = deeeeem;

%Clip all datasets to the 'active' area (where ice is present)
dem_glacier.Z(glacier_mask.Z~=1)=NaN;

%Create slope
dem_slope = gradient8(dem_glacier,'deg');
dem_slope.Z(dem_slope.Z<slope_minimum_threhsold)=slope_minimum_threhsold;

%Load x and y velocity components
ub = GRIDobj('data/velocity.tif');
surface_velocity = ub.Z.*glacier_mask.Z;

%basal slid
basal_slid = 0.25*surface_velocity;

%% Step 3: Calculate an initialize ice thickness and shape factor
Gantayat14_ice_thickness = Gantayat14_initial_icethickness(surface_velocity,basal_slid,Ac,f,p_i,g,dem_slope,n,slope_minimum_threhsold); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ice_thickness_tif = dem_glacier;        ice_thickness_tif.Z = Gantayat14_ice_thickness;       ice_thickness_tif = moving_average(ice_thickness_tif);

info = geotiffinfo([path_to_file,'dem.tif']);

geotiffwrite([result_path,'Gantayat14_ice_thickness.tif'],ice_thickness_tif.Z,deeeeem.georef.SpatialRef,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);

%% Step 4: Calculate slip law parameters
[shape_factor,flowline_xy] = cal_shape_factor(glacier_mask,outline,Gantayat14_ice_thickness);

basal_topography = deeeeem.Z-Gantayat14_ice_thickness;

dem_basal = dem_glacier;
dem_basal.Z = basal_topography;

%%% sliding parameter
[As_slid_parameter,r_basal_roughtness] = cal_sliding_parameter_As(dem_basal,B);

%%% basal shear stress
basal_shear_stress = shape_factor*p_i*g*sind(dem_slope.Z).*Gantayat14_ice_thickness; %%% Pa

%%% C
basal_slope = gradient8(dem_basal,'deg');
C = cc*nanmax(nanmax(sind(basal_slope.Z)));

%%% effective pressure
k_effect = cal_effect_k(dem_basal,shape_factor,dem_glacier,glacier_mask);
effective_pressure = k_effect*p_i*g.*Gantayat14_ice_thickness;

%%% basal sliding
basal_slid = cal_basal_slid(As_slid_parameter,C,effective_pressure,basal_shear_stress,n);
basal_slid(basal_slid<0)=0;
ice_deformation = surface_velocity-basal_slid;
ice_deformation(ice_deformation<0)=0;
mask_file = ice_deformation./ice_deformation;
basal_slid = basal_slid.*mask_file;

bu_msk = (double(glacier_mask.Z==1)+double(isnan(basal_slid)))==2;
bu_basal_slid = bu_msk*0.25.*surface_velocity;

basal_slid(isnan(basal_slid)==1) = 0;  bu_basal_slid(isnan(bu_basal_slid)==1) = 0;
basal_slid = basal_slid+bu_basal_slid;
%% Step 5: calculate ice thickness
ice_thickness = Gantayat14_initial_icethickness(surface_velocity,basal_slid,Ac,shape_factor,p_i,g,dem_slope,n,slope_minimum_threhsold);

diff_ice_thickness = nanmean(nanmean(ice_thickness - Gantayat14_ice_thickness));
i = 1;
while abs(diff_ice_thickness)>0.1
    i = i+1;
    disp(['number of iterations:',num2str(i)]);
    Gantayat14_ice_thickness = ice_thickness;

    %Calculate an initialize ice thickness and shape factor
    [shape_factor,~] = cal_shape_factor(glacier_mask,outline,Gantayat14_ice_thickness);

    %Calculate basal topography
    basal_topography = dem_glacier.Z-Gantayat14_ice_thickness;
    dem_basal = dem_glacier;
    dem_basal.Z = basal_topography;

    %%% C
    basal_slope = gradient8(dem_basal,'deg');
    C = cc*nanmax(nanmax(sind(basal_slope.Z)));

    %%% sliding parameter
    [As_slid_parameter,r_basal_roughtness] = cal_sliding_parameter_As(dem_basal,B);

    %%% basal shear stress
    basal_shear_stress = shape_factor*p_i*g*sind(dem_slope.Z).*Gantayat14_ice_thickness;

    %%% effective pressure
    k_effect = cal_effect_k(dem_basal,shape_factor,dem_glacier,glacier_mask);
    effective_pressure = k_effect*p_i*g.*Gantayat14_ice_thickness; %%% Pa

    %%% basal sliding
    basal_slid_tmp = cal_basal_slid(As_slid_parameter,C,effective_pressure,basal_shear_stress,n);
    basal_slid_tmp(basal_slid_tmp<0)=0;
    ice_deformation = surface_velocity-basal_slid_tmp;
    ice_deformation(ice_deformation<0)=0;
    mask_file = ice_deformation./ice_deformation;
    basal_slid_tmp = basal_slid_tmp.*mask_file;
    bu_msk = (double(glacier_mask.Z==1)+double(isnan(basal_slid_tmp)))==2;
    bu_basal_slid = bu_msk*0.25.*surface_velocity;
    basal_slid_tmp(isnan(basal_slid_tmp)==1) = 0;  bu_basal_slid(isnan(bu_basal_slid)==1) = 0;
    basal_slid_tmp = basal_slid_tmp+bu_basal_slid;

    basal_slid = basal_slid_tmp;

    ice_thickness = Gantayat14_initial_icethickness(surface_velocity,basal_slid,Ac,shape_factor,p_i,g,dem_slope,n,slope_minimum_threhsold);

    diff_ice_thickness = nanmean(nanmean(ice_thickness - Gantayat14_ice_thickness));

    if i>100
        break
    end

end

%% smooth result
ice_thickness_tif = dem_glacier;        ice_thickness_tif.Z = ice_thickness;       ice_thickness_tif = moving_average(ice_thickness_tif);
basal_slid_tif = dem_glacier;           basal_slid_tif.Z = basal_slid;      basal_slid_tif = moving_average(basal_slid_tif);
ice_deformation_tif = dem_glacier;      ice_deformation_tif.Z = ice_deformation;   ice_deformation_tif = moving_average(ice_deformation_tif);

ratio_ub_us = nanmean(nanmean(basal_slid./surface_velocity));

geotiffwrite([result_path,'itmsl_ice_thickness.tif'],ice_thickness_tif.Z,deeeeem.georef.SpatialRef,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite([result_path,'itmsl_basal_slid.tif'],basal_slid_tif.Z,deeeeem.georef.SpatialRef,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite([result_path,'itmsl_ice_deformation.tif'],ice_deformation_tif.Z,deeeeem.georef.SpatialRef,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);

save([result_path,'itmsl_parameters.mat'],'effective_pressure');
save([result_path,'itmsl_parameters.mat'],'As_slid_parameter','-append');
save([result_path,'itmsl_parameters.mat'],'C','-append');
save([result_path,'itmsl_parameters.mat'],'Ac','-append');
save([result_path,'itmsl_parameters.mat'],'shape_factor','-append');
save([result_path,'itmsl_parameters.mat'],'p_i','-append');
save([result_path,'itmsl_parameters.mat'],'g','-append');
save([result_path,'itmsl_parameters.mat'],'n','-append');
save([result_path,'itmsl_parameters.mat'],'flowline_xy','-append');
save([result_path,'itmsl_parameters.mat'],'r_basal_roughtness','-append');
save([result_path,'itmsl_parameters.mat'],'diff_ice_thickness','-append');
save([result_path,'itmsl_parameters.mat'],'basal_slid','-append');
save([result_path,'itmsl_parameters.mat'],'ratio_ub_us','-append');
save([result_path,'itmsl_parameters.mat'],'B','-append');
save([result_path,'itmsl_parameters.mat'],'cc','-append');
save([result_path,'itmsl_parameters.mat'],'k_effect','-append');


