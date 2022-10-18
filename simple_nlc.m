function [c, ceq] = nlc(U_vec, params, scenario)
    time = linspace(0, params("T"), params("n"));
    Fu = griddedInterpolant(time,U_vec);
    U = @(t) Fu(t); 
    
    [~, XV] = ode15s(@(t,xv) [xv(2); U(t) - params("k0") - params("k1")*xv(2) - params("k2")*xv(2)^2], time, [0; 0]);
    V = XV(:, 2); 
end 