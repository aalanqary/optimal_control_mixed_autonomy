% Specify problem params
auxdata.T = 1; % time interval in seconds
auxdata.g = 9.81; 
auxdata.k0 = auxdata.g*0.02; %k0 from the papers
auxdata.k1 = auxdata.g*1e-5;  
auxdata.k2 = auxdata.g*1e-4;
auxdata.k3  = auxdata.g*2e-4;

% Specify problem size 
auxdata.N = 10; % timesteps
auxdata.N_state = 100; 
auxdata.h = auxdata.T/auxdata.N;
auxdata.tau = linspace(0, auxdata.T, auxdata.N); % an array that splits time into spots for the control

% Specify constraints params
auxdata.eps = 0.1;
auxdata.gamma = 0.05;

% Initial guess
z = auxdata.k0 * ones(1, auxdata.N); %column vector

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
  
[U_star, info] = ipopt_auxdata(z,funcs,option);
[time_v, v] = system_solve(U_star, auxdata);

