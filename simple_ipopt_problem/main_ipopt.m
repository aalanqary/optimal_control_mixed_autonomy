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
%% Run Optimizer 
options = optimoptions('fmincon','Display','iter-detailed', ...
                        'SpecifyObjectiveGradient', true ,...
                        'SpecifyConstraintGradient', false, ...
                        'FunValCheck','on', 'DerivativeCheck', 'off',...
                        'maxfunevals',1e6, 'StepTolerance',1e-12, ...
                        'algorithm', 'interior-point');

obj = @(U) objective_gradient(U, auxdata);
funcs.objective = obj(1);
funcs.gradient = obj(2);
const = @(U) constraint_gradient(U, auxdata);
funcs.const = const(1);
% Nonlinear constraints: accepts a vector or array x and returns two arrays, c(x) and ceq(x)
funcs.jacobian = const(2);
A = [];
b = [];
% [A, b] = lc(params, scenario); 
Aeq = []; beq = []; 
[U_star,~,~,~,~,grad,~] = fmincon(fun, z, [], [], [], [], [],[],nonlcon,options);
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