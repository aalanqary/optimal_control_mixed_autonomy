function [A, b] = lc_v(auxdata)
    % Linearizatin of the constraints using ODE3 scheme
    v_0 = auxdata.v0(1);
    ts = length(auxdata.utime);

    Tv = auxdata.udt * ones(ts-1,ts);
    Tv = tril(Tv);
    Tv = [zeros(1,ts);Tv];
    %V = v_0 + Tv*U;

    A = -Tv;
    b = v_0 * ones(ts, 1);
end
    