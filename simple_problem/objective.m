function [f, df] = objective(U, auxdata)  
    [time_XV, ~, V] = system_solve(U, auxdata);
    phi_0 = 0;
    int_L_0 = -diff(V)./diff(time_XV);
    f = phi_0 + sum(int_L_0); 
    
    time = linspace(0, auxdata.T, auxdata.N);
    V = griddedInterpolant(time_XV,V, "previous"); %Todo: maybe we need to fixe the time step for forward and backward system 
    P0 = [0; 0];
    [~, P] = ode15s(@(t,p) [0; -1 + p(1) - p(2)*auxdata.k1 - p(2) * auxdata.k2*V(t)], flip(time), P0);
    df = - flip(P(:, 2));
end