function f = objective_penalty(U, auxdata)  
    [time_v, v] = system_solve(U, auxdata);
    %v = griddedInterpolant(time_v, v, "previous");

    % Extracting rho from helper functions
    pi_1 = p_1(U, v, auxdata);
    pi_2 = p_2(U, v, auxdata);

    f = -trapz(time_v, v) - auxdata.gamma * (trapz(time_v, pi_1) + trapz(time_v, pi_2));
    class(f)
end
