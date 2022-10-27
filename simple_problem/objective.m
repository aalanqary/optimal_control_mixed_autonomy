function [f, df] = objective(U, auxdata)  
    % Objective 
    tau = linspace(0, auxdata.T, auxdata.N + 1);
    [time_XV, V] = system_solve(U, auxdata);
    V = griddedInterpolant(time_XV, V, "cubic");
    phi_0 = 0;
    int_L_0 = arrayfun(@(a,b) integral(@(t) -V(t), a, b), tau(1:end-1), tau(2:end));  
    f = phi_0 + sum(int_L_0); 
    
    %Gradient
    %TODO: maybe we need to fix the time step in two system solves
    lambda0 = 0;
    [~, Lambda] = ode15s(@(t,lambda) -1 + lambda *(auxdata.k1 - 2*auxdata.k2*V(t)), flip(time_XV), lambda0);    
    Lambda = griddedInterpolant(time_XV, flip(Lambda), "cubic");
    df = arrayfun(@(a,b) integral(@(t) -Lambda(t), a, b), tau(1:end-1), tau(2:end));  
end