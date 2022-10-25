function f = objective_gradient(U, auxdata)  
    time = linspace(0, auxdata.T, auxdata.N);
    [time_XV, ~, V] = system_solve(U, auxdata);
    V = griddedInterpolant(time_XV,V, "previous"); %Todo: maybe we need to fixe the time step for forward and backward system 
    P0 = [0; 0];
    [~, P] = ode15s(@(t,p) [0; -1 + p(1) - p(2)*auxdata.k1 - p(2) * auxdata.k2*V(t)], flip(time), P0);
    f = - flip(P(:, 2));
end
