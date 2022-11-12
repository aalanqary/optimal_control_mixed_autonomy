function [time_v, v] = system_solve(U, auxdata)
    U = griddedInterpolant(auxdata.tau(1:end-1), U(1:end-1), "previous");
    time_v = linspace(0, auxdata.T, auxdata.N_state);
    v = ode5(@(t,v) U(t) - auxdata.k0 - auxdata.k1*v - auxdata.k2*v^2, time_v, 0);
end 