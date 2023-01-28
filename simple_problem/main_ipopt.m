% Specify problem params
auxdata.T = 1; % time interval in seconds
auxdata.g = 9.81; 
auxdata.k0 = auxdata.g*0.02; %k0 from the papers
auxdata.k1 = auxdata.g*1e-5;  
auxdata.k2 = auxdata.g*1e-4;
auxdata.k3  = auxdata.g*2e-4;

% Specify problem size 
auxdata.N = 10; % timesteps
auxdata.N_state = 100; % state is solving for the velocity 100 times in each of the steps for each of the time changes for the u values
auxdata.h = auxdata.T/auxdata.N;
auxdata.tau = linspace(0, auxdata.T, auxdata.N); % an array that splits time into spots for the control

% Specify constraints params
auxdata.eps = 0.1;

auxdata.gamma = 0.1;

% Initial guess - could even try starting with a one really close to the
% actual u value then keep on moving from the actual solution just
% construct a function that if less than time cutoff then high value and
% more than time then the low value
z = auxdata.k0 * ones(1, auxdata.N); %column vector
opt = optimal(1, -1, 5, 10);
% Specify functions callbakcs
funcs.objective = @objective;
funcs.gradient = @obj_grad;
funcs.constraints = @const;
% MxN matrix and needs to be sparse
% is this only finding the gradient with respect to U
funcs.jacobian = @jacobian;
% contains 1s at values which Jacobian can have nonzero entries
funcs.jacobianstructure = @jacobianstructure;

% don't need these - bounds for the variables
option.lb = -50*ones(size(z)) ; % Lower bound on the variables.
option.ub = 50*ones(size(z)) ; % Upper bound on the variables.

% The constraint functions are bounded to zero = 0;
option.cl =  [0, -auxdata.gamma, -auxdata.gamma];
option.cu = [0, inf, inf];

% Set the IPOPT options -These might have to change
option.ipopt.tol                   = 1e-1;
option.ipopt.max_iter             = 5;
option.ipopt.print_level           = 5;
%option.ipopt.jac_c_constant        = 'yes'; % indicates whether all
%equality constraints are linear also jac_d_constant for inequality
%constraints
option.ipopt.hessian_approximation = 'limited-memory'; % limited-memory quasi-Newton approximation
option.ipopt.mu_strategy           = 'adaptive'; % update strategy for barrier parameter

% Set up the auxiliary data.
option.auxdata = auxdata;

% Specify problem params
aux.T = 1;
aux.g = 9.81;
aux.k0 = aux.g*0.02;
aux.k1 = aux.g*1e-5;  
aux.k2 = aux.g*1e-4;
aux.k3  = aux.g*2e-4;

% Specify problem size 
aux.N = 10;
aux.N_state = 100; 
aux.h = aux.T/aux.N;
aux.tau = linspace(0, aux.T, aux.N);
% Specify constraints params
aux.eps = 0.5;
aux.gamma = 0.5;
nonlcon = @(U) constraint_gradient(U, aux);
options = optimoptions('fmincon','Display','iter-detailed', ...
                        'SpecifyObjectiveGradient', true ,...
                        'SpecifyConstraintGradient', false, ...
                        'FunValCheck','on', 'DerivativeCheck', 'off',...
                        'maxfunevals',1e6, 'StepTolerance',1e-12, ...
                        'algorithm', 'interior-point');

fun = @(U) objective_gradient(U, aux);
A = [];
b = [];
% [A, b] = lc(params, scenario); 
Aeq = []; beq = []; 
[U_star,~,~,~,~,grad,~] = fmincon(fun, opt, [], [], [], [], [],[],nonlcon,options);
[time_v, v] = system_solve(U_star, aux);
%[U_star, info] = ipopt_auxdata(opt,funcs,option); % just change this to fmincon and add in their parameters
%[time_v, v] = system_solve(U_star, auxdata);

function [f, df] = objective_gradient(U, aux)
    f = objective(U, aux); 
    df = obj_grad(U, aux);
end 

function [c, ceq, dc, dceq] = constraint_gradient(U, auxdata)
    [c, ceq] = constfmincon(U, auxdata); 
    if nargout > 2 
        [dc, dceq] = const_grad(U, auxdata);
    end 
    
end 
