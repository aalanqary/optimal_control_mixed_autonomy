function j = J(X, V, U, scenario, params) % both X and U are matrices
%     VU = [V, U];
    VU = V;
    idx = scenario("I_a");
    X_diff = [];
    des_v = params("des_v");
    if idx(1) == 1
        xl = scenario("x_leader");
        X_diff = xl(params("t_int")) - X(:, 1) - params("l") - params("safe_dist");
        idx = idx(2:end);
    end 
    
    X_diff = [X_diff, X(:, idx-1) - X(:, idx) - params("l") - params("safe_dist")];
    
    integrand = sum((VU-des_v(params("t_int"))).^2, 2) + sum(log(min(1, max(X_diff, eps))).^2, 2);
    %     j = diff(params("t_int"))'* integrand(1:end-1);
    j = params("T")/params("nt") * sum(integrand)
    if min(X_diff) <= 0
        display(params("T")/params("nt") * sum(sum((VU-des_v(params("t_int"))).^2, 2)), 'j_velocity')
        display(params("T")/params("nt") * sum(sum(log(min(1, max(X_diff, eps))).^2, 2)), 'j_position')
%         display(j, 'j_position')
    end 
end 