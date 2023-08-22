function [A, b] = lc(auxdata)
    Xl = auxdata.xl(auxdata.utime);
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

    A = [(Tx*Tv + auxdata.d_min*Tv); -(Tx*Tv + auxdata.d_max*Tv);  -Tv];
    b = [Xl - (x_0 + v_0 * sum(Tx, 2)) - auxdata.l - auxdata.d_min*v_0;
        - Xl + (x_0 + v_0 * sum(Tx, 2)) + auxdata.l + auxdata.d_max*v_0;
        v_0 * ones(ts, 1)];
end
    