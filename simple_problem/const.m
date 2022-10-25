function [c, ceq] = const(U, auxdata)  
    [time_XV, ~, V] = system_solve(U, auxdata);
    ceq = V(end);
    c = [U - auxdata.g - auxdata.k3 * V.^2;
        - U - auxdata.g - auxdata.k3 * V.^2];

%     [time_XV, ~, V] = system_solve(U, auxdata);
%     U_t = ; %todo: Make the vector U have the same dim as the vector V
%     
%     phi_1 = (V(end)).^2; 
%     h1 = auxdata.g + auxdata.k3*V.^2 - U_t; % >= 0 
%     h1 = h1 * (h1 < auxdata.eps) - ((h - auxdata.eps).^2/3*auxdata.eps) * ((-auxdata.eps <= h) .* ())
% 
%     h2 =  U_t - auxdata.g - auxdata.k3*V.^2; % >= 0 
%     h = h * (h < auxdata.eps) - ((h - auxdata.eps).^2/3*eps) * ((-auxdata.eps <= h) .* ())
%     int_L_1 = diff() ./ diff(time_XV); 
%     f1 = phi_1 + sum(int_L+1);
end