function [X, V, A] = system_solve_b(b, vl, xl, time, x0) 
%%IMPORTANT NOTE: this only solve systems that have AV in index 1 then all
%%HV and time spacing of AV is 1 s
    xlf = griddedInterpolant(time, xl);
    vlf = griddedInterpolant(time, vl);
    F = @(t, XV) [XV(2); ACC_b(b, xlf(t), XV(1), vlf(t), XV(2))];

    XV = ode5(F, time, [x0; vl(1)]);
    V = XV(:, 2);
    X = XV(:, 1);

    if nargout > 2
        A = ACC_b(b, xl, X, vl, V); 
    end 
end

function a = ACC_b(b, xl, xf, vl, vf)
    alpha = b(1); 
    beta = b(2); 
    k = b(3); 
    safe_dist = b(4); 

    l = 5; 
    v_max = 35;

    head_way = xl - xf - l;
    c = tanh(l + safe_dist);
    Vopt = v_max .* ((tanh(k*head_way - safe_dist) + c)./(1 + c));
    a = alpha *(Vopt - vf) + ...
        beta * (vl - vf) ./ (head_way.^2);
    a;
end 
