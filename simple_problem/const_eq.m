function [c, ceq] = const_eq(U, auxdata)  
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

    % Using Traditional Approach
    if true
        % First constraint (eq): v(T) = 0  
        ceq(1) = v(end); 
    end 

end