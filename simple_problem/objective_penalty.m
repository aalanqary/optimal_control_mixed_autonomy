function f = objective_penalty(U, auxdata)  
    [time_v, v] = system_solve(U, auxdata);

%     figure(1)
%     plot(time_v, v)
%     drawnow
%     figure(2)
%     plot(auxdata.tau, U)
%     plot(U)
%     drawnow

    disp(U)
    disp("Original Constraint Violation: ")
    disp(const(U, auxdata));
    [h1, h2] = const(U, auxdata);
    disp(h1);
    disp(h2);

    %Plot original constraints
    %plot(auxdata.tau, h1, 'g');
    %drawnow;

%     integ_h1 = sum(min(h1, 0).^2);
%     integ_h2 = sum(min(h2, 0).^2);

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
    p1 = p_1(U_1, v, auxdata);
    % Plotting/testing smooth penalty approximation
    figure(1);
    plot(time_v, p1, 'g', auxdata.tau, h1, 'r');
    legend('ρ1', 'h1')
    title("ρ1 vs h1")
    xlabel("t")
    ylabel("Violation")
    drawnow;


    
    % Extracting rho from helper functions
    p2 = p_2(U_1, v, auxdata);
    % Plotting/testing smooth penalty approximation
    figure(2);
    plot(time_v, p2, 'g', auxdata.tau, h2, 'r');
    legend('ρ2', 'h2')
    title("ρ2 vs h2")
    xlabel("t")
    ylabel("Violation")
    drawnow;

    obj = -trapz(time_v, v); 
    constraint_penalty = (trapz(time_v, p1) + trapz(time_v, p2));

    % Display values
    disp("Original Objective Value: " + obj)
    disp("Constraint Value: " + constraint_penalty)
    disp("Total Constraint Penalty: " + (-auxdata.gamma * constraint_penalty))

    f = obj - auxdata.gamma * constraint_penalty;
    disp("Total Objective Value: " + f)
    %trapz(time_v, v)
    %trapz(time_v, p1)
    %trapz(time_v, p2)
end
