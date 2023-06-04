function [dc, diff_vec] = finite_diff(U, h, auxdata)
    [X, V, A, Fx, Fv, Fa, Fu] = constraint_gradient_min_computation(U, auxdata);
    [utimelength, unum] = size(U);
    dc_all = zeros(unum*utimelength, 2*unum);
    c_all = zeros(unum * 2, 1); 
    for i = 1:length(auxdata.Ia)
        [ci, dci] = constraint_gradient_min_max_index(X, V, A, Fx, Fv, Fa, Fu, auxdata, i, U);
        c_all(2*i-1:2*i) = ci;
        dc_all(:, 2*i-1:2*i) = dci;
    end
    c = c_all;
    dc = dc_all;
    diff_vec = zeros(unum*utimelength, 2*unum);
    for i = 1:utimelength
        for m = 1:unum
            U_shift = U;
            U_shift(i, m) = U(m, 1) + h;
            [X, V, A, Fx, Fv, Fa, Fu] = constraint_gradient_min_computation(U_shift, auxdata);
            [timelength, num] = size(V);
            c_all = zeros(unum * 2, 1); 
            for j = 1:length(auxdata.Ia)
                [ci] = constraint_min_max_index_vector(X, V, A, Fx, Fv, Fa, Fu, auxdata, j, U_shift);
                c_all(2*j-1:2*j) = ci;
            end
            c_shift = c_all;
            % these with respect to the first AV
            for k = 1:2*unum
                diff_vec(i + (m - 1) * utimelength, k) = (c(k) - c_shift(k)) / h; % should be for 2 * unum values
            end
        end
    end

%     total_diff = zeros(2*unum, 1);
%     for i = 1:unum
%         min_vec_1 = diff_vec(1:utimelength, 2 * i - 1);
%         min_vec_2 = diff_vec(1:utimelength, 2 * i);
%         max_vec_1 = diff_vec(utimelength + 1:2*utimelength, 2 * i - 1);
%         max_vec_2 = diff_vec(utimelength + 1:2*utimelength, 2 * i);
%         total_diff(2 * i - 1) = trapz(auxdata.utime, min_vec_1) + trapz(auxdata.utime, min_vec_2);
%         total_diff(2 * i) = trapz(auxdata.utime, max_vec_1) + trapz(auxdata.utime, max_vec_2);
%     end
        
% %     diff_vec = (c - c_shift) / h;
%     finite_diff = zeros(4, 1);
%     grad_summation = zeros(4, 1);
%     total_diff = zeros(4, 1);
%     for i = 1:2
%         min_vec = diff_vec(1:timelength, i);
%         max_vec = diff_vec(timelength + 1:2*timelength, i);
%         finite_diff(2 * i - 1) = trapz(auxdata.time, min_vec);
%         finite_diff(2 * i) = trapz(auxdata.time, max_vec);
%         min_grad = dc(1:utimelength, i);
%         max_grad = dc(utimelength + 1:2*utimelength, i);
%         grad_summation(2 * i - 1) = trapz(auxdata.utime, min_grad);
%         grad_summation(2 * i) = trapz(auxdata.utime, max_grad);
%         total_diff(2 * i - 1) = grad_summation(2 * i - 1) - finite_diff(2 * i - 1);
%         total_diff(2 * i) = grad_summation(2 * i) - finite_diff(2 * i);
%     end
       





