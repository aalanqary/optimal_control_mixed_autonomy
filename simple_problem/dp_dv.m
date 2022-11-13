function [dp1_dv,dp2_dv] = dp_dv(u, v, auxdata)
%D_RHO_D_V Returns
dp1_dv = f_adj_2(u, v, auxdata);
dp2_dv = f_adj_3(u, v, auxdata);
end

function [f] = f_adj_2(u, v, auxdata)
    h = auxdata.g + auxdata.k3*v^2 - u; 
    if h < -auxdata.eps
        f = 2*auxdata.k3*v;
    elseif (h >= -auxdata.eps) && (h<= auxdata.eps)
        f = (-1/4*auxdata.eps) * (4*auxdata.k3^2 * v^3 + 4*auxdata.k3*v*(auxdata.g - u - auxdata.eps));
    else
        f = 0; 
    end
end 

function [f] = f_adj_3(u, v, auxdata)
    h = auxdata.g + auxdata.k3*v^2 + u; 
    if h < -auxdata.eps
        f = 2*auxdata.k3*v(t);
    elseif (h >= -auxdata.eps) && (h<= auxdata.eps)
        f = (-1/4*auxdata.eps) * (4*auxdata.k3^2 * v^3 + 4*auxdata.k3*v*(auxdata.g + u - auxdata.eps));
    else
        f = 0; 
    end
end 

