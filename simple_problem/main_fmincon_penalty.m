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
auxdata.gamma = 1;
auxdata.gamma_const = 0.2;

z = auxdata.k0 * ones(1, auxdata.N); %column vector
%z = sin(auxdata.tau); Not required

%% Run Optimizaer 
options = optimoptions('fmincon','Display','iter-detailed', ...
                        'SpecifyObjectiveGradient', true ,...
                        'SpecifyConstraintGradient', false, ...
                        'FunValCheck','on', 'DerivativeCheck', 'off',...
                        'maxfunevals',1e6, 'StepTolerance',1e-12, ...
                        'algorithm', 'interior-point');

fun = @(U) objective_penalty_gradient(U, auxdata);
% Nonlinear constraints: accepts a vector or array x and returns two arrays, c(x) and ceq(x)
nonlcon = @(U) constraint_gradient(U, auxdata);
A = [];
b = [];
A_min = -100*ones(10);
A_max = 100*ones(10);
% [A, b] = lc(params, scenario); 
Aeq = []; beq = []; 
[U_star,~,~,~,~,grad,~] = fmincon(fun, z, [], [], [], [], A_min,A_max,nonlcon,options);
[time_v, v] = system_solve(U_star, auxdata);

% Example U violation test
U = [9.9455, 8.9526, 10.5842, 5.6858, 0.4102, -4.8671, -9.2073, -9.8769, -9.5680, 0.1962];
const(U, auxdata)


function [f, df] = objective_penalty_gradient(U, auxdata)
    f = objective_penalty(U, auxdata); 
    df = obj_grad_penalty(U, auxdata);
end 

function [c, ceq] = constraint_gradient(U, auxdata)
    [c, ceq] = const_eq(U, auxdata); 
%     if nargout > 
%         [dc, dceq] = const_eq_grad(U, auxdata);
%     end 
    
end 