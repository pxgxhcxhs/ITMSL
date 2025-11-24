function [ice_thickness] = Gantayat14_initial_icethickness(surface_velocity,basal_slid,A0,f,p_i,g,dem_slope,n,slope_minimum_threhsold)
%Calculate icethickness

%Relevant citations to read are:

%Gantayat, Prateek, Anil V. Kulkarni, and J. Srinivasan. "Estimation of
%Ice Thickness Using Surface Velocities and Slope: Case Study at Gangotri Glacier, India."
%Journal of Glaciology 60, no. 220 (2014): 277-82. https://dx.doi.org/10.3189/2014JoG13J078.

%Do not allow slopes lower than minimum threshold
dem_slope.Z(dem_slope.Z<slope_minimum_threhsold)=slope_minimum_threhsold; 

%velocity:m/s; A:Pa^-3*s^-1; ice_density: kg/m^3,g: m/s^2,h: m
ice_thickness = (((n+1)*(surface_velocity-basal_slid)/(60*60*24*365))./...
    (2*A0*(f*p_i*g*sind(dem_slope.Z)).^n)).^(1/(n+1));

end

