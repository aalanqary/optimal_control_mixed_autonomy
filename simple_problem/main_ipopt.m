% Specify problem params
auxdata.T = 1;
auxdata.g = 9.81;
auxdata.k0 = auxdata.g*0.02;
auxdata.k1 = auxdata.g*1e-5;  
auxdata.k2 = auxdata.g*1e-4;
auxdata.k3  = auxdata.g*2e-4;

% Specify problem size 
auxdata.N = 10;
auxdata.N_state = 100; 
auxdata.h = auxdata.T/auxdata.N;
auxdata.tau = linspace(0, auxdata.T, auxdata.N);

% Specify constraints params
auxdata.eps = 0;
auxdata.gamma = 0;

% Initial guess
z = auxdata.k0 * ones(1, auxdata.N); %column vector

% Specify functions callbakcs
funcs.objective = @objective;
funcs.gradient = @obj_grad;
funcs.constraints = @const;
% MxN matrix and needs to be sparse
% is this only finding the gradient with respect to U
funcs.jacobian = @jacobian;
% 10x10 matrix need to change this so that the 7 columns are zeros
%figure out how to set m and n
funcs.jacobianstructure = @jacobianstructure;

% don't need these - bounds for the variables
option.lb = -Inf*ones(size(z)) ; % Lower bound on the variables.
option.ub = Inf*ones(size(z)) ; % Upper bound on the variables.

% The constraint functions are bounded to zero = 0;
option.cl =  [0, -inf, -inf];
option.cu = [0, 0, 0];

% Set the IPOPT options -These might have to change
option.ipopt.print_level           = 3;
%option.ipopt.jac_c_constant        = 'yes';
option.ipopt.hessian_approximation = 'limited-memory';
option.ipopt.mu_strategy           = 'adaptive';
option.ipopt.tol                   = 1e-7;
% Set up the auxiliary data.
option.auxdata = auxdata;
  
[U_star, info] = ipopt_auxdata(z,funcs,option);
[time_v, v] = system_solve(U_star, auxdata);

