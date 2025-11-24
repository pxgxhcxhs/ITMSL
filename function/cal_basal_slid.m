function [basal_slid] = cal_basal_slid(As_slid_parameter,C,effective_pressure,basal_shear_stress,n)
%Relevant citations to read are:

%O., Gagliardini, D., Cohen, P., R?back, and, T., and Zwinger. 
%"Finite-Element Modeling of Subglacial Cavities and Related Friction Law." 
%Journal of Geophysical Research Earth Surface  (2007).

%Schoof, C. "The Effect of Cavitation on Glacier Sliding." 
%Proceedings of the Royal Society A: Mathematical, Physical and Engineering Science 461, no. 2055 (2005): 609-27.

%Helanow, Christian, Neal R. Iverson, Jacob B. Woodard, and Lucas K. Zoet. 
%"A Slip Law for Hard-Bedded Glaciers Derived from Observed Bed Topography." 
%Science Advances 7, no. 20 (2021): eabe7798. https://dx.doi.org/doi:10.1126/sciadv.abe7798.

%Gimbert, F., A. Gilbert, O. Gagliardini, C. Vincent, and L. Moreau. 
%"Do Existing Theories Explain Seasonal to Multi©\Decadal Changes in Glacier Basal Sliding Speed?", 
%Geophysical Research Letters 48, no. 15 (2021). https://dx.doi.org/10.1029/2021gl092858.

% Conversion units
As_slid_parameter = As_slid_parameter/10^(6*n); %m/(yr*Pa^n)

%Calculate basal sliding velocity Ub(m/yr):
basal_slid = As_slid_parameter*(C*effective_pressure.*basal_shear_stress).^n./(C^n*effective_pressure.^n-basal_shear_stress.^n);

end

