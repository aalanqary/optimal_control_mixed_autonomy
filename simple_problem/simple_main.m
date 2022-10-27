% Specify problem params
auxdata.T = 10;
auxdata.g = 9.81;
auxdata.k0 = auxdata.g*0.02;
auxdata.k1 = auxdata.g*1e-5;
auxdata.k2 = auxdata.g*1e-4;
auxdata.k3  = auxdata.g*2e-4;

% Specify problem size 
auxdata.N = 100;
auxdata.h = auxdata.T/auxdata.N;

% Specify IPOPT options
% options.auxdata = auxdata;


z = 0 * ones(auxdata.N, 1); %column vector
% z = [9.81 * ones(auxdata.N/2, 1); -9.81 * ones(auxdata.N/2, 1)];

%% Run Optimizaer 

options = optimoptions('fmincon','Display','iter-detailed', ...
                        'SpecifyObjectiveGradient', true,...
                        'FunValCheck','on', 'DerivativeCheck', 'off',...
                        'maxfunevals',1e6, 'StepTolerance',1e-12, ...
                        'algorithm', 'interior-point');

fun = @(U) objective(U, auxdata);
% Nonlinear constraints: accepts a vector or array x and returns two arrays, c(x) and ceq(x)
nonlcon = [];%@(U) const(U, auxdata);
A = [];
b = [];
% [A, b] = lc(params, scenario); 
Aeq = []; beq = []; 
lb = []; ub = [];
[U_star,~,~,~,~,grad,~] = fmincon(fun, z, A, b, Aeq, beq, lb,ub,nonlcon,options);
% a_min = -3 * ones(size(U0)); 
% a_max = 1.5 * ones(size(U0));
% 
% funcs.objective = @(U) U;
% funcs.gradient = fun;
% funcs.constraints = @(nonlcol) nonlcol;
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

[x, info] = ipopt(x0,funcs,option);
%[U_star, f_val, ~, output, ~, grad] = fmincon(fun, U0, A, b, Aeq, beq, a_min, a_max, nonlcon, options);

% [U_star,~,~,~,~,grad,~] = fmincon(fun, U0, A,b,Aeq,beq,lb,ub,nonlcon,options);
%info;

% these need to be uncommented and added in later
%Fu = griddedInterpolant(scenario("time"),U_star);
%U_star = @(t) Fu(t);
%[X_star, V_star] = system_solve(U_star, params, scenario);

