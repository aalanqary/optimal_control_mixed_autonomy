function f = objective_penalty(U, auxdata)  
    [time_v, v] = system_solve(U, auxdata); 
    U = griddedInterpolant(auxdata.tau, U, "previous");
    U_1 = U(time_v);
    %v = v(auxdata.tau);
    if isrow(v)
        v = v';
    end 
    if isrow(U_1)
        U_1 = U_1'; 
    end 

    % Extracting rho from helper functions
    pi_1 = p_1(U_1, v, auxdata);
    %size(pi_1)
    pi_2 = p_2(U_1, v, auxdata);

    f = -trapz(time_v, v); 
    f = f - auxdata.gamma * (trapz(time_v, pi_1) + trapz(time_v, pi_2));
    %trapz(time_v, v)
    %trapz(time_v, pi_1)
    %trapz(time_v, pi_2)
end
