function XV_dot = F(t, XV, U, scenario, params) 
    X = XV(1:length(scenario("config")));
    V = XV(length(scenario("config")) + 1: end);
    X_dot = V;
    V_dot = zeros(length(V), 1);
    
    xl = scenario("x_leader");
    xl = xl(t);
    vl = scenario("v_leader");
    vl = vl(t);
    X_l = [xl; X]; 
    V_l = [vl; V]; 
    a = ACC(X_l(scenario("I_h")), X(scenario("I_h")), ... 
                                 V_l(scenario("I_h")), V(scenario("I_h")), ...
                                 params);
    V_dot(scenario("I_h")) = a;
    V_dot(scenario("I_a")) = U; 
    XV_dot = [X_dot; V_dot];
end 

