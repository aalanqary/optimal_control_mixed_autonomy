%% Define problem 
platoon = [1, 0, 1, 0];
data_dir = 'data_v2_preprocessed_west';
data_file = dir(sprintf('%s/*7050.csv', data_dir));

scenario = prepare_data(data_dir, data_file.name, 2, platoon, 1000);

%specify parameters
params = containers.Map;
params("ode_opt") = odeset('RelTol',1e-8,'AbsTol',1e-10); % []
% IDM params
params("safe_dist") = 2.5; %float
params("v_max") = 35; %float
params("l") = 5; %float
params("alpha") = 0.5;
params("beta") = 21;
% objective function params 
params("mu1") = 1;
params("mu2") = 1;
% optimizer params 
params("eps") = 2;
params("gamma") = 120;
params("use_gradient") = true; 
params("initalize") = "IDM"; 

%Initial solution  
[X0, V0, U0] = initial_solution(containers.Map(keys(scenario), values(scenario)), containers.Map(keys(params), values(params)));

%% Run Optimizaer 

%options = optimoptions('fmincon','Display','iter-detailed', ...
                        %'SpecifyObjectiveGradient',params("use_gradient"),...
                        %'FunValCheck','on', 'DerivativeCheck', 'off',...
                        %'maxfunevals',1e6, 'StepTolerance',1e-12, ...
                        %'algorithm', 'interior-point');

fun = @(U) objective_gradient(U, params, scenario); 
% Nonlinear constraints: accepts a vector or array x and returns two arrays, c(x) and ceq(x)
nonlcon = @(U) nlc(U, params, scenario);
A = [];
b = [];
% [A, b] = lc(params, scenario); 
Aeq = []; beq = []; 
a_min = -3 * ones(size(U0)); 
a_max = 1.5 * ones(size(U0));

funcs.objective = @(U) U;
funcs.gradient = fun;
funcs.constraints = @(nonlcol) nonlcol;
% Jacobian
%funcs.jacobian  = @(nonlcol) jacobian(nonlcol, U);
%funcs.jacobian = @train_NLP_constraints_jacobian;
%funcs.jacobianstructure = @train_NLP_constraints_jacobian_pattern;
%funcs.jacobianstructure = @(nonlcol) jacobian(nonlcol, U);
%starting point
x0 = X0;


% Lower/Upper bounds on constraints (don't know if these are correct
option.cl = a_min;    
option.cu = a_max;

% Set the IPOPT options -These might have to change
option.ipopt.print_level           = 3;
option.ipopt.jac_c_constant        = 'yes';
option.ipopt.hessian_approximation = 'limited-memory';
option.ipopt.mu_strategy           = 'adaptive';
option.ipopt.tol                   = 1e-7;

[x info] = ipopt(x0,funcs,option);
%[U_star, f_val, ~, output, ~, grad] = fmincon(fun, U0, A, b, Aeq, beq, a_min, a_max, nonlcon, options);

% [U_star,~,~,~,~,grad,~] = fmincon(fun, U0, A,b,Aeq,beq,lb,ub,nonlcon,options);
%info;

% these need to be uncommented and added in later
%Fu = griddedInterpolant(scenario("time"),U_star);
%U_star = @(t) Fu(t);
%[X_star, V_star] = system_solve(U_star, params, scenario);

