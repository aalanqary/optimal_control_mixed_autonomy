function [c, ceq] = const(U, auxdata)  
    [time_XV, ~, V] = system_solve(U, auxdata); % Given the control input, solve for X and V
    ceq = V(end);
    c = [U - auxdata.g - auxdata.k3 * V(1:end-1).^2;
        - U - auxdata.g - auxdata.k3 * V(1:end-1).^2];

    U_t = U; %todo: Make the vector U have the same dim as the vector V

    phi_1 = (V(end)).^2; 
    h1 = auxdata.g + auxdata.k3*V.^2 - U_t; % >= 0  would this be c[2]
    bound_h1 = (-auxdata.eps <= h1) + (auxdata.eps >= h1);
    h1 = h1 * (h1 < auxdata.eps) - ((h1 - auxdata.eps).^2/4*auxdata.eps) * (bound_h1 >= 2);

    h2 =  U_t - auxdata.g - auxdata.k3*V.^2; % >= 0 would this be c[1]
    bound_h2 = (-auxdata.eps <= h2) + (auxdata.eps >= h2);
    h2 = h2 * (h2 < auxdata.eps) - ((h2 - auxdata.eps).^2/4*auxdata.eps) * (bound_h2 >= 2);
    
    
    int_L_1 = diff(h1) ./ diff(time_XV); 
    f1 = phi_1 + sum(int_L_1+1);

    int_L_2 = diff(h2) ./ diff(time_XV); 
    f2 = phi_1 + sum(int_L_2+1);

    % [Nonlinear inequality constraints, Nonlinear equality constraints]
    % intial c = [U - auxdata.g - auxdata.k3 * V.^2;
       % - U - auxdata.g - auxdata.k3 * V.^2];
    c(1) = f1;
    c(2) = f2;
end