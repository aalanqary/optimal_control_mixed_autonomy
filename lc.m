function [A, b] = lc(auxdata)
    Xl = auxdata.xl(auxdata.utime);
    x_0 = auxdata.x0(1);
    v_0 = auxdata.v0(1);
    ts = length(auxdata.utime);
    T = auxdata.utime(end)/ts *  ones(ts-1,ts);
    T = tril(T);
    T = [zeros(1,ts);T];

    %V = v_0 * ones(ts, 1) + T*U;
    %X = x_0 * ones(ts, 1) + v_0 * auxdata.utime + T*T*U;
    A = [T*T; -T*T;  -T];
    b = [Xl - (x_0 * ones(ts, 1) + v_0 * auxdata.utime) - auxdata.l - auxdata.d_min;
        - Xl + (x_0 * ones(ts, 1) + v_0 * auxdata.utime) + auxdata.l + auxdata.d_max;
        v_0 * ones(ts, 1)];
end
    