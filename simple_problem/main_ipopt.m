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

z = auxdata.k0 * ones(1, auxdata.N); %column vector
z = sin(auxdata.tau);
a_min = -3 * ones(size(z)); 
a_max = 1.5 * ones(size(z));
b_min = -3 * ones(size(z)); 
b_max = 1.5 * ones(size(z));
%% Run Optimizer 

funcs.objective = @(U) U;
funcs.gradient = @(U) obj_grad(U, auxdata);
%
funcs.const = @(U) const(U, auxdata);

% Nonlinear constraints: accepts a vector or array x and returns two arrays, c(x) and ceq(x)
grad = @(U) const_grad(U, auxdata);
% MxN matrix and needs to be sparse
funcs.jacobian = @(U) sparse(jacobian(U, auxdata));
display(funcs.jacobian)
display(size(funcs.jacobian))
%funcs.jacobianstructure = @() sparse(ones(size(funcs.jacobian)));
% 10x10 matrix need to change this so that the 7 columns are zeros
%figure out how to set m and n
funcs.jacobianstructure = @() sparse(ones(3, 10));

A = [];
b = [];
% [A, b] = lc(params, scenario); 
Aeq = []; beq = []; 

% Lower/Upper bounds on constraints - these might not be correct
option.cl = a_min;    
option.cu = a_max;

% Lower/Upper bounds on variables - these might not be correct
option.lb = b_min;    
option.ub = b_max;

% Set the IPOPT options -These might have to change
option.ipopt.print_level           = 3;
%option.ipopt.jac_c_constant        = 'yes';
option.ipopt.hessian_approximation = 'limited-memory';
option.ipopt.mu_strategy           = 'adaptive';
option.ipopt.tol                   = 1e-7;

[U_star, info] = ipopt_auxdata(z,funcs,option);
[time_v, v] = system_solve(U_star, auxdata);


function [f, df] = objective_gradient(U, auxdata)
    f = objective(U, auxdata); 
    df = obj_grad(U, auxdata);
end 

function [c, ceq, dc, dceq] = constraint_gradient(U, auxdata)
    [c, ceq] = const(U, auxdata); 
    if nargout > 2 
        [dc, dceq] = const_grad(U, auxdata);
    end 
    
end 