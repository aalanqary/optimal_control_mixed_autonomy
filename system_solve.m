function [X, V, A] = system_solve(U, params, scenario_simple)

    [~, XV] = ode15s(@(t,XV) F(t, XV, U(t)', scenario_simple, params), scenario_simple("time"), [scenario_simple("x0"); scenario_simple("v0")], params("ode_opt"));
    V = XV(:, end - length(scenario_simple("config")) + 1: end);
    X = XV(:, 1:end - length(scenario_simple("config")));
    
    if nargout > 2
        A = zeros(size(X));
        xl = scenario_simple("xl");
        xl = xl(scenario_simple("time"));
        vl = scenario_simple("vl");
        vl = vl(scenario_simple("time"));
        Xl = [xl, X]; 
        Vl = [vl, V]; 
        Ah = ACC(Xl(:, scenario_simple("Ih")), X(:, scenario_simple("Ih")), ... 
                      Vl(:, scenario_simple("Ih")), V(:, scenario_simple("Ih")), params); 
        A(:, scenario_simple("Ih")) = Ah;
        A(:, scenario_simple("Ia")) = U(scenario_simple("time")); 
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

function XV_dot = F(t, XV, U, scenario_simple, params) 
    X = XV(1:length(scenario_simple("config")));
    V = XV(length(scenario_simple("config")) + 1: end);
    
    xl = scenario_simple("xl");
    xl = xl(t);
    vl = scenario_simple("vl");
    vl = vl(t);
    Xl = [xl; X]; 
    Vl = [vl; V]; 
    a = ACC(Xl(scenario_simple("Ih")), X(scenario_simple("Ih")), ... 
                                 Vl(scenario_simple("Ih")), V(scenario_simple("Ih")), ...
                                 params);
    X_dot = V;
    V_dot = zeros(length(V), 1);
    V_dot(scenario_simple("Ih")) = a;
    V_dot(scenario_simple("Ia")) = U; 
    XV_dot = [X_dot; V_dot];
end 