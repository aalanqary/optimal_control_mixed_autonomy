function [X0, V0, U0] = initial_solution(scenario, params)
      if params("initalize") == "IDM"
          original_I_a = scenario("Ia");
          scenario("config") = zeros(1, length(scenario("config")));
          scenario("Ia") = find(scenario("config"));
          scenario("Ih") = find(scenario("config") - 1);
          params("l") = 10 * params("l");
          U_temp = @(t) []; 
          [X0, V0, U0] = system_solve(U_temp, params, scenario);
          U0 = U0(:, original_I_a);
      else
          error("Initialization method not found")
      end 
%       Fu = griddedInterpolant(scenario("time"),U0);
%       U0 = @(t) Fu(t);
end 