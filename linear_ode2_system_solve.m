function [X, V] = linear_ode2_system_solve(U, auxdata)
    x_0 = auxdata.x0(1);
    v_0 = auxdata.v0(1);
    ts = length(auxdata.utime);
    T = auxdata.udt * ones(ts,ts);
    T = T - diag(diag(T)./2);
    T(:, 1) = T(:, 1)./2;
    T = tril(T);
    T(1,1) = 0; 
    V = v_0 + T*U;
%     X = ones(ts, 1) + T*V; 
    X = x_0 + v_0 * sum(T,2) + T*T*U;
end 