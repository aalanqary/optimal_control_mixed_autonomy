function E = fuel_consumption(X, V, U, scenario, params)

    % Energy term 
    A = zeros(size(X));
    xl = scenario("x_leader");
    xl = xl(params("t_int"));
    vl = scenario("v_leader");
    vl = vl(params("t_int"));
    X_l = [xl, X]; 
    V_l = [vl, V]; 
    A(:, scenario("I_h")) = ACC(X_l(:, scenario("I_h")), X(:, scenario("I_h")), ... 
              V_l(:, scenario("I_h")), V(:, scenario("I_h")), params);
%     a = diff(V(:, scenario("I_h")))./ diff(params("t_int"));
%     a = [a; a(end, :)];
%     A(:, scenario("I_h")) = a;
    A(:, scenario("I_a")) = U;
    E = simplified_fuel_model(V,A,'RAV4');
    
%     % Objective function value 
%     fc = params("T")/params("nt") * sum(E, "all");
    
end 