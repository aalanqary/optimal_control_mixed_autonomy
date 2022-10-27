function [time_XV, V] = system_solve(U, auxdata)
    tau = linspace(0, auxdata.T, auxdata.N+1);
    sys_time = linspace(0, auxdata.T, 10 * auxdata.N); 
    U = griddedInterpolant(tau(1:end-1), U, "previous");
    [time_XV, XV] = ode15s(@(t,xv) U(t) - auxdata.k0 - auxdata.k1*xv - auxdata.k2*xv^2, sys_time, 0);
    V = XV;
end 