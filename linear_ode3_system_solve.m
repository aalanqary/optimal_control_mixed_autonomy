function [X, V] = linear_ode3_system_solve(U, auxdata)
    x_0 = auxdata.x0(1);
    v_0 = auxdata.v0(1);
    ts = length(auxdata.utime);
    Tv = auxdata.udt * ones(ts-1,ts);
    Tv = tril(Tv);
    Tv = [zeros(1,ts);Tv];
    V = v_0 + Tv*U;

    Tx = auxdata.udt * ones(ts,ts);
    Tx = Tx - diag(diag(Tx)./2);
    Tx(:, 1) = Tx(:, 1)./2;
    Tx = tril(Tx);
    Tx(1,1) = 0; 
    
%     X = x_0 + Tx*V; 
    X = x_0 + sum(Tx, 2) * v_0 + Tx*Tv*U; 
end 