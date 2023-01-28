function [c, ceq] = constfmincon(U, auxdata)  
    % c(x, u) <= 0 
    % ceq(x, u) = 0
    c = [];
    [time_v, v] = system_solve(U, auxdata); 
    v = griddedInterpolant(time_v, v, "previous");
    v = v(auxdata.tau);
    
    if isrow(v)
        v = v';
    end 
    if isrow(U)
        U = U'; 
    end 
    
    if false
        % First constraint (eq): v(T) = 0  
        ceq(1) = v(end); 
        
        c = [U - auxdata.g - auxdata.k3 * v.^2; % Second constraint (ineq): u(t) <= g + k3 * v(t)
            - U - auxdata.g - auxdata.k3 * v.^2]; % Third constraint (ineq): u(t) >= -g - k3 * v(t)
    end 

    if false
        % First constraint (eq): v(T) = 0  
        ceq(1) = v(end); 
        % Second constraint (ineq): u(t) <= g + k3 * v(t)
        h1 = auxdata.g + auxdata.k3 * v - U;
        ceq(2) = sum(min(h1, 0).^2);
        % Third constraint (ineq): u(t) >= -g - k3 * v(t)
        h2 =  auxdata.g + auxdata.k3*v.^2 + U; 
        ceq(3) = sum(min(h2, 0).^2);

        c = [];
    end 

    % Using Traditional Approach
    if true
        % First constraint (eq): v(T) = 0  
        c(1) = v(end);
        % Second constraint (ineq): u(t) <= g + k3 * v(t)
        % g + k3 * v(t) - u(t) >= 0
        pi_1 = p_1(U, v, auxdata);
        ceq(2) = trapz(auxdata.tau, pi_1);
        % Third constraint (ineq): u(t) >= -g - k3 * v(t)
        % u(t) + g + k3 * v(t) >= 0
        pi_2 = p_2(U, v, auxdata);
        ceq(3) = trapz(auxdata.tau, pi_2);
    end 
end