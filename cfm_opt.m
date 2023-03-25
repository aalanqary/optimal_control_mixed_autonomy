T = 704; 
auxdata.dt = 0.1;
auxdata.udt = 1; 
auxdata.utime = (0:auxdata.udt:T)';
time = (0:auxdata.dt:T)';

optdata = readtable("optimal_data/2021-04-22-12-47-13_2T3MWRFVXLW056972_masterArray_0_7050.csv");
vl = optdata.Var3;
% vl = griddedInterpolant(time,vl);
xl = optdata.Var1;
xf = optdata.Var2;
vopt = optdata.Var4;
aopt = optdata.Var5;

x0 = xl(1) - xf(1); 
% v0 = vl(1); 

options = optimoptions('fmincon','Display','iter-detailed', ...
                        'SpecifyObjectiveGradient', false ,...
                        'FunValCheck','on', 'DerivativeCheck', 'off',...
                        'maxfunevals',1e6, 'StepTolerance',1e-12, ...
                        'algorithm', 'interior-point', ...
                        'ConstraintTolerance',1e-6);

fun = @(b) objective(b, xl, vl, time, x0, vopt); 
nonlcon = []; %@(U) nlc(U, auxdata);
A = [];
b = [];
% [A, b] = lc(auxdata); 
Aeq = []; beq = []; 
a_min = [00.0001,10,0.1,2.5]; %-3 * ones(size(U0)); 
a_max = [0.1, 600, 1, 100]; %3 * ones(size(U0));
tic
[U_star, f_val, ~, output, ~, grad] = fmincon(fun, [0.1000, 525.0000, 0.2000, 4.0000], A, b, Aeq, beq, a_min, a_max, nonlcon, options);
toc


function z = objective(b, xl, vl, time, x0, aopt)
    
    % Objective 
    [X, V, A] = system_solve_b(b, vl, xl, time, x0);
    z = sum((aopt - A).^2)./time(end);

    if isnan(z)
        z = inf;
    end 
end 