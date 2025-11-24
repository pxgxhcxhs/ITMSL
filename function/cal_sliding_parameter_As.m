function [As_slid_parameter,r_basal_roughtness] = cal_sliding_parameter_As(dem_basal,B)
%Relevant citations to read are:

%O., Gagliardini, D., Cohen, P., R?back, and, T., and Zwinger. 
%"Finite-Element Modeling of Subglacial Cavities and Related Friction Law." 
%Journal of Geophysical Research Earth Surface  (2007).

%Calculate sliding parameters As:
[m,n] = size(dem_basal.Z);
extra_value = zeros(m,n);


%The search window is 3x3
for i=2:m-1
    for j=2:n-1
        tmp_array = [];
        tmp_array = dem_basal.Z(i-1:i+1,j-1:j+1);
        max_value = max(max(tmp_array));
        min_value = min(min(tmp_array));
        if tmp_array(2,2)==max_value
            extra_value(i,j) = 1;
        elseif tmp_array(2,2)==min_value
            extra_value(i,j) = 2;
        end
    end
end

m1=[];n1=[];m2=[];n2=[];
[m1,n1] = find(extra_value==1);
[m2,n2] = find(extra_value==2);

for i=1:length(m1)
    %%% calculate wavelength
    dist = [];
    dist=(sqrt((m1(i,1)-m1).^2+(n1(i,1)-n1).^2))*dem_basal.cellsize; 
    [wave_length(i,1),~]=min(dist(dist~=0));
    
    %%% calculate amplitude
    dist = [];
    dist=(sqrt((m1(i,1)-m2).^2+(n1(i,1)-n2).^2))*dem_basal.cellsize;
    [~,ind]=min(dist(dist~=0));
    amplitude(i,1) = dem_basal.Z(m1(i,1),n1(i,1))-dem_basal.Z(m2(ind,1),n2(ind,1));    
    
end

wave_length_ave = nanmean(wave_length);
amplitude_ave = nanmean(amplitude)/2;   %%%%%%%% pxg20231211

r_basal_roughtness = amplitude_ave/wave_length_ave;

As_slid_parameter = B*wave_length_ave*(-2*10^(-5)+0.013*r_basal_roughtness+0.0262*r_basal_roughtness^2)/r_basal_roughtness^2;

end

