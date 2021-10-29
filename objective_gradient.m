function [z, dz] = objective_gradient(U, scenario, params)
    
    global fun_interp
    % Solve forward system 
    U_func = @(t) fun_interp(U, t);
    options_ode = [];
    [~, XV] = ode45(@(t,XV) F(t, XV, U_func(t)', scenario, params), params("t_int"), [scenario("x_0"); scenario("v_0")], options_ode);
    
    V = XV(:, end - length(scenario("I_h")) + 1: end);
    X = XV(:, 1:end - length(scenario("I_h")));
    
    % Evaluate Objective 
    z = J(X, V, U_func(params("t_int")), scenario, params);
    
    if 1==1
        figure(10)
        x_leader = scenario("x_leader");
        plot(params("t_int"),x_leader(params("t_int")),params("t_int"),X)
        drawnow
        
        figure(11)
        v_leader = scenario("v_leader");
        plot(params("t_int"),v_leader(params("t_int")),params("t_int"),V)
        drawnow
    end    
%         
    
    if nargout >= 2
        display("solving backward")
        % Solve backward system
        [~, dim] = size(XV);
        P_0 = zeros(dim, 1);
        [~,PQ] = ode45(@(t,PQ) F_adjoint(t, PQ, X, V, U, scenario, params), flip(params("t_int")), P_0, options_ode);

        PQ = flip(PQ,1);
        Q = PQ(:, end - length(scenario("I_h")) + 1: end);
        P = PQ(:, 1:end - length(scenario("I_h")));

        dz = [];
        for i = scenario("I_a")
            grad = P(:, i) - L_partial(params("t_int"), i, X, V, U, "v_a", scenario, params);
            if ismember(i+1, scenario("I_h"))
                grad = grad + ACC_partial(X(:, i), X(:,i+1), U(:, scenario("I_a")==i), V(:, scenario("I_h")==i+1), 3, params) .* Q(:, scenario("I_h")==i+1); 
            end
            dz = [dz, -1*grad];
        end 
        display("Done solving backward")
    end 
%     if 1==1
%         figure(2)
%         plot(params("t_int"),grad)
%         drawnow
%     end   

end 