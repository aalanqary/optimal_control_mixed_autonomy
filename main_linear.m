platoon = [1];
[auxdata, leader] = problem_auxdata(platoon, "linear", "simple");

%Initial solution - smooth leader trajectory
U0 = diff(smoothdata(leader.v(auxdata.utime), 'movmean', 10)) ./ diff(auxdata.utime);
U0 = [U0;0];
[X0, V0, A0] = system_solve(U0, auxdata, leader);

%% Run Optimizaer 

options = optimoptions('fmincon','Display','iter-detailed', ...
                        'SpecifyObjectiveGradient', true ,...
                        'SpecifyConstraintGradient', false ,...
                        'FunValCheck','on', 'DerivativeCheck', 'off',...
                        'maxfunevals',1e6, 'StepTolerance',1e-12, ...
                        'algorithm', 'sqp');

fun = @(U) objective_gradient_accel(U, auxdata, leader); 
nonlcon = []; 
[A, b] = lc_time_headway(auxdata, leader); 
Aeq = []; beq = []; 
a_min = [];
a_max = []; 
[U_star, f_val, ~, output, ~, grad] = fmincon(fun, U0, A, b, Aeq, beq, a_min, a_max, nonlcon, options);
[X_star, V_star, A_star] = system_solve(U_star, auxdata, leader);

