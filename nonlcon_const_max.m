function [c, ceq, dc, dceq] = nonlcon_const_max(U, auxdata)
    [X, V, A, Fx, Fv, Fa, Fu] = constraint_gradient_min_computation(U, auxdata);
    [timelength, num] = size(U);
    c_all = zeros(num, 1*2);
    dc_all = zeros(num*timelength, num*2);
    for i = 1:length(auxdata.Ia)
        [ci, dci] = constraint_gradient_min_max_index(X, V, A, Fx, Fv, Fa, Fu, auxdata, i, U);
        c_all(:, i) = ci;
        dc_all(:, 2*i-1:2*i) = dci; % address this as a sparse matrix
    end
    ceq = [];
    dceq = [];
    c = c_all;
    dc = dc_all;
    figure(1);
    plot(c);
    figure(2);
    plot(dc);
end