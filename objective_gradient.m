function [z, dz] = objective_gradient(U, params, scenario)
    
    % Objective 
    [X, V] = system_solve(U, params, scenario);
    z = J(X, V, U, scenario, params);
    
    % Gradient
    if nargout >= 2
        % Solve backward system
        display("solving backward")
        [~, Q] = back_system_solve(X, V, U, params, scenario);
        display("Done solving backward")
        % Compute gradeint
        dz = -1* Q(:, scenario("I_a")) +  L_partial(X(:, scenario("I_a")), V(:, scenario("I_a")), U, "u", params);
    end 
    
    
    if 1==2
        figure(10)
        x_leader = scenario("x_leader");
        plot(params("t_int"),x_leader(params("t_int")),params("t_int"),X)
        ylabel("position")
        drawnow
        
        figure(11)
        v_leader = scenario("v_leader");
        plot(params("t_int"),v_leader(params("t_int")),params("t_int"), V)
        ylabel("velocity")
        drawnow
    end    
end 