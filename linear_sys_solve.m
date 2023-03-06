function [X, V] = linear_sys_solve(U, auxdata)
    x_0 = auxdata.x0(1);
    v_0 = auxdata.v0(1);
    ts = length(auxdata.utime);
    T = (auxdata.utime(end)+1)/ts *  ones(ts-1,ts);
    T = tril(T);
    T = [zeros(1,ts);T];
    
    V = v_0 * ones(ts, 1) + T*U;
%     X = x_0 * ones(ts, 1) + T*V; 
    X = x_0 * ones(ts, 1) + v_0 * auxdata.utime + T*T*U;
end 