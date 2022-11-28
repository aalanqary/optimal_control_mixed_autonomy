function [p_2] = p_2(U, v, auxdata)
    h2 =  auxdata.g + auxdata.k3*v.^2 + U; 
    bound_h2 = (-auxdata.eps <= h2) + (auxdata.eps >= h2);
    p_2 = h2 .* (h2 < -auxdata.eps) - ((h2 - auxdata.eps).^2/4*auxdata.eps) .* (bound_h2 >= 2);
end

