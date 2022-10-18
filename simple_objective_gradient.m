% page 157
function [z, dz] = objective_gradient(U_vec, params, scenario)
    
    % Objective 
    time = linspace(0, params("T"), params("n"));
    Fu = griddedInterpolant(time,U_vec);
    U = @(t) Fu(t); 
    
    [~, XV] = ode15s(@(t,xv) [xv(2); U(t) - params("k0") - params("k1")*xv(2) - params("k2")*xv(2)^2], time, [0; 0]);
    V = XV(:, 2); 
    z = params("T")/params("n") * sum(V, "all"); 
    
    % Gradient
    if nargout >= 2
        %Fx = griddedInterpolant(scenario("time"),X);
        Fv = griddedInterpolant(scenario("time"),V);
        [~, P] = ode15s(@(t,p) [0; -1 + p(1) - p(2)*k1 - p(2) * k2*Fv(t)], time, [0; 0]);
        dz = -P(:, 2);
    end 
 
end
