function [c, ceq, dc, dceq] = nonlcon_const(U, auxdata)
    [X, V, A, Fx, Fv, Fa, Fu, PQ0, P_dot, Q_dot] = constraint_gradient_min_computation(U, auxdata);
    [timelength, num] = size(U);
    c_all = zeros(num, 1);
    dc_all = zeros(num, 2 * length(auxdata.utime));
    for i = 1:num
        [ci, dci] = constraint_gradient_min_index(X, V, A, Fx, Fv, Fa, Fu, PQ0, P_dot, Q_dot, auxdata, i, U(:, i));
        c_all(i) = ci;
        dc_all(i, :) = [dci', dci']; % address this as a sparse matrix
    end
    ceq = [];
    dceq = [];
    c = c_all;
    dc = dc_all';
end