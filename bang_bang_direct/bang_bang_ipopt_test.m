%% Compute NLP solution with IPOPT
%
% This function computes the numerical solution of the optimal control problem 
% as a direct transcription problem (NLP) and unsing IPOPT Optimizer.
% 
% See the paper for model details.
%
% Input parameters:
%  N                number of grid points
%  p_data           structure with problem data
%    p_data.g       acceleration of gravity g      ;
%    p_data.T_size  time horizon
%    p_data.k0      constant friction (normalized with mass)
%    p_data.k1      friction linearly dependent on speed (normalized with mass)
%    p_data.k2      drag friction (normalized with mass)
%    p_data.k3      down force (normalized with mass)
%  compute_accuracy compute Norm-2 w.r.t. exact solution: (true/false)
%                   Only for k2 = 0, k3 = 0
% plotting          plot solution: (true/false)
%
% Return:
%  norm2:          is empty if  compute_accuracy=false

% No formal inputs for now
function [z,elapsed,ok] = bang_bang_ipopt_test()

  % 10 samples for test
  N = auxdata.N

  % Total time
  T = 20

  lb = [ -Inf*ones(2*(N+1),1) ; -Inf*ones(N,1) ] ; % Lower bound constraints to -infinity
  ub = [  Inf*ones(2*(N+1),1) ;  Inf*ones(N,1) ] ;
  cl = [ zeros(2*N+3,1) ; zeros(2*N,1) ] ; % 2*N+3 equality constraints
  cu = [ zeros(2*N+3,1) ; Inf*ones(2*N,1) ] ;

  options.lb = lb ; % Lower bound on the variables.
  options.ub = ub ; % Upper bound on the variables.

  % The constraint functions are bounded to zero
  options.cl = cl ;
  options.cu = cu ;

  % Set up the auxiliary data.
  options.auxdata = auxdata ;
  
  % Set the IPOPT options.
  options.ipopt.jac_d_constant   = 'no';
  options.ipopt.hessian_constant = 'yes';
  options.ipopt.mu_strategy      = 'adaptive';
  options.ipopt.max_iter         = 400;
  options.ipopt.tol              = 1e-10;
  
  % The callback functions.
  funcs.objective         = @objective;
  funcs.gradient          = @gradient;
  funcs.constraints       = @constraints;
  funcs.jacobian          = @direct_method_constraints_jacobian;
  funcs.jacobianstructure = @direct_method_constraints_jacobian_pattern;
  if use_hessian
    funcs.hessian           = @direct_method_hessian;
    funcs.hessianstructure  = @direct_method_hessian_pattern;
  else
    options.ipopt.hessian_approximation      = 'limited-memory';
    options.ipopt.limited_memory_update_type = 'bfgs' ; % {bfgs}, sr1 = 6; % {6}
    options.ipopt.limited_memory_max_history = 50 ;
  end
  %options.ipopt.derivative_test = 'first-order';
  %options.ipopt.derivative_test = 'second-order';
  
  % Initial guess solution -> Set discretization at each time to 0
  x = zeros(N+1,1) ;
  v = zeros(N+1,1) ;
  u = zeros(N,1) ;
  z0 = [x ;v ;u ] ;

  tic
  [z, info] = ipopt_auxdata(z0,funcs,options);
  elapsed = toc ;

  ok = info.iter ;

end
