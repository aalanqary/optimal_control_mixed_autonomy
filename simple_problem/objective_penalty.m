function f = objective_penalty(U, auxdata)  
    [time_v, v] = system_solve(U, auxdata);

    % Extracting rho from helper functions
    p_1 = rho_1()
    p_2 = rho_2()

    f = -trapz(time_v, v) - auxdata.gamma * (trapz(auxdata.tau, p_1) + trapz(auxdata.tau, p_2));
end
