function f = objective_penalty(U, auxdata)  
    [time_v, v] = system_solve(U, auxdata);

    % Recomputation of pi_1
    h1 = auxdata.g + auxdata.k3 * v - U;
    bound_h1 = (-auxdata.eps <= h1) + (auxdata.eps >= h1);
    pi_1 = h1 .* (h1 < -auxdata.eps) - (((h1 - auxdata.eps).^2)/4*auxdata.eps) .* (bound_h1 >= 2);

    % Recomputation of pi_2
    h2 =  auxdata.g + auxdata.k3*v.^2 + U; 
    bound_h2 = (-auxdata.eps <= h2) + (auxdata.eps >= h2);
    pi_2 = h2 .* (h2 < -auxdata.eps) - ((h2 - auxdata.eps).^2/4*auxdata.eps) .* (bound_h2 >= 2);

    f = -trapz(time_v, v) - auxdata.gamma * (trapz(auxdata.tau, phi_1) + trapz(auxdata.tau, phi_2));
end
