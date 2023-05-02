function [c, ceq, dc, dceq] = nonlcon_const(U, auxdata)
    [X, V, A, Fx, Fv, Fa, Fu] = constraint_gradient_min_computation(U, auxdata);
    [timelength, num] = size(U);
    c_all = zeros(num, 1);
    dc_all = zeros(num*length(auxdata.utime), num);
    for i = 1:length(auxdata.Ia)
        [ci, dci] = constraint_gradient_min_index(X, V, A, Fx, Fv, Fa, Fu, auxdata, i, U);
        c_all(:, i) = ci;
        dc_all((i-1)*length(auxdata.utime)+1:i*length(auxdata.utime), :) = dci; % address this as a sparse matrix
    end
    ceq = [];
    dceq = [];
    c = c_all;
    dc = dc_all;
end