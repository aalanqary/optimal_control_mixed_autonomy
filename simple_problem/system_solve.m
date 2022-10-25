function [time_XV, X, V] = system_solve(U, auxdata)
    time = linspace(0, auxdata.T, auxdata.N);
    U = griddedInterpolant(time, U, "previous");
    [time_XV, XV] = ode15s(@(t,xv) [xv(2); U(t) - auxdata.k0 - auxdata.k1*xv(2) - auxdata.k2*xv(2)^2], time, [0; 0]);
    V = XV(:, 2); 
    X = XV(:, 1); 
end 