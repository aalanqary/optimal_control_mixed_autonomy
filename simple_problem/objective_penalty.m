function f = objective_penalty(U, auxdata)  
    [time_v, v] = system_solve(U, auxdata);

%     figure(1)
%     plot(time_v, v)
%     drawnow
%     figure(2)
%     plot(auxdata.tau, U)
%     plot(U)
%     drawnow

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

    obj = -trapz(time_v, v); 
    constraint_penalty = (trapz(time_v, pi_1) + trapz(time_v, pi_2));

    % Display values
    disp("Original Objective Value: " + obj)
    disp("Constraint Value: " + constraint_penalty)
    disp("Total Constraint Penalty: " + (-auxdata.gamma * constraint_penalty))

    f = obj - auxdata.gamma * constraint_penalty;
    disp("Total Objective Value: " + f)
    %trapz(time_v, v)
    %trapz(time_v, pi_1)
    %trapz(time_v, pi_2)
end
