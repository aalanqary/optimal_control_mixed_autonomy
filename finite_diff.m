function [c_shift, c, finite_diff, total_diff] = finite_diff(U, h, auxdata)
    [X, V, A, Fx, Fv, Fa, Fu] = constraint_gradient_min_computation(U, auxdata);
    [vtimelength, num] = size(U);
    dc_all = zeros(2*vtimelength, num*num);
    [timelength, num] = size(V);
    c_all = zeros(2*timelength, length(auxdata.Ia));
    for i = 1:length(auxdata.Ia)
        [ci, dci] = constraint_gradient_min_max_index(X, V, A, Fx, Fv, Fa, Fu, auxdata, i, U);
        c_all(2*i-1:2*i) = ci;
        dc_all(:, 2*i-1:2*i) = dci;
    end
    c = c_all;
    dc = dc_all;

    U_shift = U + h;
    [X, V, A, Fx, Fv, Fa, Fu] = constraint_gradient_min_computation(U_shift, auxdata);
    [timelength, num] = size(V);
    c_all = zeros(2*timelength, length(auxdata.Ia));
    for i = 1:length(auxdata.Ia)
        [ci] = constraint_gradient_min_max_index_vector(X, V, A, Fx, Fv, Fa, Fu, auxdata, i, U_shift);
        c_all(:, i) = ci;
    end
    c_shift = c_all;
    
    diff_vec = (c - c_shift) / h;
    finite_diff = zeros(4, 1);
    grad_summation = zeros(4, 1);
    total_diff = zeros(4, 1);
    for i = 1:2
        min_vec = diff_vec(1:timelength, i);
        max_vec = diff_vec(timelength + 1:2*timelength, i);
        finite_diff(2 * i - 1) = trapz(auxdata.time, min_vec);
        finite_diff(2 * i) = trapz(auxdata.time, max_vec);
        min_grad = dc(1:vtimelength, i);
        max_grad = dc(vtimelength + 1:2*vtimelength, i);
        grad_summation(2 * i - 1) = trapz(auxdata.utime, min_grad);
        grad_summation(2 * i) = trapz(auxdata.utime, max_grad);
        total_diff(2 * i - 1) = grad_summation(2 * i - 1) - finite_diff(2 * i - 1);
        total_diff(2 * i) = grad_summation(2 * i) - finite_diff(2 * i);
    end
       





