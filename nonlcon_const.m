function [c, ceq, dc, dceq] = nonlcon_const(U, auxdata)
    [X, V, A, Fx, Fv, Fa, Fu, PQ0, P_dot, Q_dot] = constraint_gradient_min_computation(U, auxdata);
    c_all = zeros(length(auxdata.time), length(U));
    dc_all = zeros(length(auxdata.time), length(U));
    for i = 1:length(U_vec)
        [c, dc] = constraint_gradient_min_index(X, V, A, Fx, Fv, Fa, Fu, PQ0, P_dot, Q_dot, auxdata, i);
        c_all(:, i) = c;
        dc_all(:, i) = dc;
    end
    ceq = [];
    dceq = [];
    c = c_all;
    dc = dc_all;
end