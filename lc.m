function [A, b] = lc(params, scenario)
    x_0 = scenario("x_0");
    x_0 = x_0(1); 
    v_0 = scenario("v_0");
    v_0 = v_0(1); 
    X_l = scenario("x_leader");
    X_l = X_l(params("t_int"));
    T = (params("T")/params("nt")) * ones(params("nt")-1,params("nt"));
    T = [zeros(1,params("nt"));T];
    T = tril(T);
%     V = v_0 * ones(params("nt"), 1) + T*U;
%     X = x_0 * ones(params("nt"), 1) + v_0 * params("t_int") + T*T*U;
    A = [T*T; -T]; 
    b = [X_l - x_0 * ones(params("nt"), 1) - v_0 * params("t_int") - params("eps")* ones(params("nt"), 1); ...
        v_0 * params("t_int")];
end