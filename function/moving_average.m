function [slope_sm] = moving_average(slope)
%Smoothing the data
[m,n] = size(slope.Z);
slope_sm = slope;
for s=1:5
    for i=2:m-1
        for j=2:n-1
            mid_point = slope_sm.Z(i,j);
            a = numel(mid_point(isnan(mid_point)));
            
            round_point(1,1) = slope_sm.Z(i-1,j);
            round_point(1,2) = slope_sm.Z(i+1,j);
            round_point(1,3) = slope_sm.Z(i,j-1);
            round_point(1,4) = slope_sm.Z(i,j+1);
            b = numel(round_point(isnan(round_point)));
            if a==1 || b==4
                                
            else
                slope_sm.Z(i,j) = 0.5*slope_sm.Z(i,j)+0.125*4/(4-b)*(nansum(round_point));
            end
            
        end
    end
end
end

