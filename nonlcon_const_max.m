function [c, ceq, dc, dceq] = nonlcon_const_max(U, auxdata)
    [X, V, A, Fx, Fv, Fa, Fu] = constraint_gradient_min_computation(U, auxdata);
    [timelength, num] = size(U);
    c_all = zeros(num * 2, 1); % should be a vector of 2 * num in format of car 1 min, car 1 max, car 2 min...
    dc_all = zeros(2*timelength, num*num); % this is basically for number of cars * utime for each constraint
    for i = 1:length(auxdata.Ia)
        [ci, dci] = constraint_gradient_min_max_index(X, V, A, Fx, Fv, Fa, Fu, auxdata, i, U);
        c_all(2*i-1:2*i) = ci;
        dc_all(:, 2*i-1:2*i) = dci; % address this as a sparse matrix
    end
    ceq = [];
    dceq = [];
    c = c_all;
    dc = dc_all;
%     display(c);
%       figure(1);
%       title("Constraint Gradient")
%       plot(dc);
%       legend('car 1 min', 'car 2 min', 'car 1 max', 'car 2 max');
%       title("Constraint Gradients")
%      drawnow;
end