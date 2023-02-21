%% Define problem 
global fun_eval
fun_eval = 0;
auxdata.platoon = [1, 0];
auxdata.len_platoon = length(auxdata.platoon);
auxdata.Ia = find(auxdata.platoon);
auxdata.Ih = find(auxdata.platoon - 1);
% data_dir = 'data_v2_preprocessed_west';
% data_file = dir(sprintf('%s/*7050.csv', data_dir));

% Leader's trajectory
auxdata.utime = (0:1:50)';
auxdata.time = (0:0.1:50)';
auxdata.vl = @(t) 30 + 0.*t; 
x0 = 0.27 * auxdata.vl(0) * flip(0:1:auxdata.len_platoon)'; 
auxdata.v0 = ones(auxdata.len_platoon, 1) * auxdata.vl(0); 
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[~, xl] = ode45(@(t,x) auxdata.vl(t), auxdata.time, x0(1), opts);
auxdata.x0 = x0(2:end);
auxdata.xl = griddedInterpolant(auxdata.time, xl);

% Bando-FtL params
auxdata.safe_dist = 2.5; 
auxdata.v_max = 35; 
auxdata.alpha = 0.5;
auxdata.beta = 21;
auxdata.l = 5; 

% objective function auxdata 
auxdata.mu1 = 1;
auxdata.mu2 = 1;

% Constraints auxdata
auxdata.d_min = auxdata.safe_dist;

% optimizer auxdata 
auxdata.eps = 2;
auxdata.gamma = 120;

%Initial solution  
% [X0, V0, U0] = initial_solution(containers.Map(keys(scenario), values(scenario)), containers.Map(keys(auxdata), values(auxdata)));
U0 = zeros(length(auxdata.utime), length(auxdata.Ia));
[X0, V0] = system_solve(U0, auxdata);

%% Run Optimizaer 

options = optimoptions('fmincon','Display','iter-detailed', ...
                        'SpecifyObjectiveGradient',true ,...
                        'FunValCheck','on', 'DerivativeCheck', 'off',...
                        'maxfunevals',1e6, 'StepTolerance',1e-12, ...
                        'algorithm', 'interior-point');

fun = @(U) objective_gradient_vel(U, auxdata); 
nonlcon = @(U) nlc(U, auxdata);
A = [];
b = [];
% [A, b] = lc(auxdata, scenario); 
Aeq = []; beq = []; 
a_min = -3 * ones(size(U0)); 
a_max = 3 * ones(size(U0));
[U_star, f_val, ~, output, ~, grad] = fmincon(fun, U0, A, b, Aeq, beq, a_min, a_max, nonlcon, options);
[X_star, V_star] = system_solve(U_star, auxdata);

