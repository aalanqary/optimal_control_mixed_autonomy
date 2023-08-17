function [A, b] = lc_time_headway(auxdata, leader)
    % Linearizatin of the constraints using ODE3 scheme
    Xl = leader.x(auxdata.utime);
    x_0 = auxdata.x0(1);
    v_0 = auxdata.v0(1);
    ts = length(auxdata.utime);

    Tv = auxdata.udt * ones(ts-1,ts);
    Tv = tril(Tv);
    Tv = [zeros(1,ts);Tv];
    %V = v_0 + Tv*U;

    Tx = auxdata.udt * ones(ts,ts);
    Tx = Tx - diag(diag(Tx)./2);
    Tx(:, 1) = Tx(:, 1)./2;
    Tx = tril(Tx);
    Tx(1,1) = 0; 
    %X = x_0 + sum(Tx, 2) * v_0 + Tx*Tv*U; 

    A = [(Tx*Tv + auxdata.h_min*Tv); -(Tx*Tv + auxdata.h_max*Tv);  -Tv];
    b = [Xl - (x_0 + v_0 * sum(Tx, 2)) - auxdata.l - auxdata.h_min*v_0;
        - Xl + (x_0 + v_0 * sum(Tx, 2)) + auxdata.l + auxdata.h_max*v_0;
        v_0 * ones(ts, 1)];
end
    