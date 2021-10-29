function j = J(X, V, U, scenario, params) % both X and U are matrices
%     VU = [V, U];
    load(['RAV4_coeffs.mat'])       % load vehicle coefficients
    idx = scenario("I_a");
    X_diff = [];
    des_v = params("des_v");
    if idx(1) == 1
        xl = scenario("x_leader");
        X_diff = xl(params("t_int")) - X(:, 1) - params("l") - params("safe_dist");
        idx = idx(2:end);
    end 
    A = [];
    for i = scenario("I_h")
           x = X(:, i);
           v = V(:, scenario("I_h")==i);
        if i == 1
           xl = scenario("x_leader");
           xl = xl(params("t_int"));
           vl = scenario("v_leader");
           vl = vl(params("t_int"));
        elseif ismember(i-1, scenario("I_h"))
           xl = X(:, i-1);
           vl = V(:, scenario("I_h")==i-1);
        else 
           xl = X(:, i-1);
           vl = U(:, scenario("I_a")==i-1);
        end 
        A = [A, ACC(xl, x, vl, v, params)];
    end 
    if size(V) ~= size(A)
        error("Velocity and acceleration sizes are not the same")
    end 
    X_diff = [X_diff, X(:, idx-1) - X(:, idx) - params("l") - params("safe_dist")];
    j_velocity = sum((V-des_v(params("t_int"))).^2, 2);
    j_penalize = sum(log(min(1, max(X_diff, eps))).^2, 2);
    j_energy = sum(simplified_fuel_model(V,A,'RAV4'), 2);
    j = params("T")/params("nt") * sum(j_velocity + j_penalize + j_energy);
%     if min(X_diff) <= 0
%         display(params("T")/params("nt") * sum(sum((VU-des_v(params("t_int"))).^2, 2)), 'j_velocity')
%         display(params("T")/params("nt") * sum(sum(log(min(1, max(X_diff, eps))).^2, 2)), 'j_position')
% %         display(j, 'j_position')
%     end 
end 