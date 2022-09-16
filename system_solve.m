function [X, V, A] = system_solve(U, params, scenario)

    [~, XV] = ode15s(@(t,XV) F(t, XV, U(t)', scenario, params), scenario("time"), [scenario("x0"); scenario("v0")], params("ode_opt"));
    V = XV(:, end - length(scenario("config")) + 1: end);
    X = XV(:, 1:end - length(scenario("config")));
    
    if nargout > 2
        A = zeros(size(X));
        xl = scenario("xl");
        xl = xl(scenario("time"));
        vl = scenario("vl");
        vl = vl(scenario("time"));
        Xl = [xl, X]; 
        Vl = [vl, V]; 
        Ah = ACC(Xl(:, scenario("Ih")), X(:, scenario("Ih")), ... 
                      Vl(:, scenario("Ih")), V(:, scenario("Ih")), params); 
        A(:, scenario("Ih")) = Ah;
        A(:, scenario("Ia")) = U(scenario("time")); 
    end 
    
%     Linear solver
%     x_0 = scenario("x_0");
%     x_0 = x_0(1); 
%     v_0 = scenario("v_0");
%     v_0 = v_0(1); 
%     T = (params("T")/params("nt")) * ones(params("nt")-1,params("nt"));
%     T = [zeros(1,params("nt"));T];
%     T = tril(T);
%     V = v_0 * ones(params("nt"), 1) + T*U;
%     X = x_0 * ones(params("nt"), 1) + v_0 * params("t_int") + T*T*U;

end 

function XV_dot = F(t, XV, U, scenario, params) 
    X = XV(1:length(scenario("config")));
    V = XV(length(scenario("config")) + 1: end);
    
    xl = scenario("xl");
    xl = xl(t);
    vl = scenario("vl");
    vl = vl(t);
    Xl = [xl; X]; 
    Vl = [vl; V]; 
    a = ACC(Xl(scenario("Ih")), X(scenario("Ih")), ... 
                                 Vl(scenario("Ih")), V(scenario("Ih")), ...
                                 params);
    X_dot = V;
    V_dot = zeros(length(V), 1);
    V_dot(scenario("Ih")) = a;
    V_dot(scenario("Ia")) = U; 
    XV_dot = [X_dot; V_dot];
end 