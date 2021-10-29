function XV_dot = F(t, XV, U, scenario, params) %X \in R(2M+N) and U \in R(N) 
    V = XV(end - length(scenario("I_h")) + 1: end);
    X = XV(1:end - length(scenario("I_h")));
    
    X_dot = zeros(length(X), 1);
    X_dot(scenario("I_h")) = V; 
    X_dot(scenario("I_a")) = U; 
    
    V_dot = zeros(length(V), 1);
    for i = scenario("I_h")
        x = X(i);
        v = V(scenario("I_h") == i);
        if i == 1
            xl = scenario("x_leader");
            xl = xl(t);
            vl = scenario("v_leader");
            vl = vl(t);
        elseif ismember(i - 1, scenario("I_h"))
            vl = V(scenario("I_h") == i-1);
            xl = X(i-1);
        else
            vl = U(scenario("I_a") == i-1);
            xl = X(i-1);
        end
        
        V_dot(scenario("I_h")==i) = ACC(xl, x, vl, v, params);
    end     
    XV_dot = [X_dot; V_dot];
end 

