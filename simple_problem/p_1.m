function [p_1] = p_1(U, v, auxdata)
    %size(v)
    %size(U)
    h1 = auxdata.g + auxdata.k3 * v.^2 - U;
    bound_h1 = (-auxdata.eps <= h1) + (auxdata.eps >= h1);
    p_1 = h1 .* (h1 < -auxdata.eps) - (((h1 - auxdata.eps).^2)/4*auxdata.eps) .* (bound_h1 >= 2);
end