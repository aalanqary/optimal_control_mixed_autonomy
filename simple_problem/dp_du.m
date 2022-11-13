function [dp1_du,dp2_du] = dp_du(U, v, auxdata)
%DP_DU Summary of this function goes here
%   Detailed explanation goes here
dp1_du = f_int_2(u, v, auxdata);
dp2_du = f_int_3(u, v, auxdata);
end

function f = f_int_2(u, v, auxdata)
    h = auxdata.g + auxdata.k3*v.^2 - u; 
    bound_h = (-auxdata.eps <= h) + (auxdata.eps >= h);
    
    f = -1 * (h < -auxdata.eps) + ...
        (-(1/4*auxdata.eps) * (2 * u - 2*(auxdata.g + auxdata.k3*v.^2 - auxdata.eps))) .* (bound_h >= 2); 
end 


function f = f_int_3(u, v, auxdata)
    h = auxdata.g + auxdata.k3*v.^2 + u; 
    bound_h = (-auxdata.eps <= h) + (auxdata.eps >= h);
    f = 1 * (h < -auxdata.eps) + ...
        (-(1/4*auxdata.eps) * (2 * u + 2*(auxdata.g + auxdata.k3*v.^2 - auxdata.eps))) .* (bound_h >= 2); 
end 

