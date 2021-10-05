%% Define the scenario
% first element is leading car, 0 represents human driver and 1 autonomous
scenario = containers.Map;
scenario("config") = [1, 0, 0,1];
scenario("I_a") = find(scenario("config"));
scenario("I_h") = find(scenario("config") - 1);
scenario("v_leader") = @(t) 2;
scenario("x_leader") = @(t) 4*t + 60;
x_leader = scenario("x_leader");
scenario("x_0") = 2*[20; 15;10; 5];
scenario("v_0") = [0; 0];

params = containers.Map;
params("t_0") = 0;
params("T") = 10;
params("nt") = 200;
params("t_int") = linspace(params("t_0"), params("T"), params("nt"))'; 
params("u_0") = @(t) [0*t + 3, 0*t + 3];
params("safe_dist") = 2; %float
params("v_max") = 40; %float
params("l") = 3; %float
params("alpha") = 0.1;
params("beta") = 10;
params("des_v") = 10; %float

useGradient = true;

global fun_interp coefs
fun_interp = @(val, t) interp1(params("t_int"), val, t);
coefs = containers.Map;
coefs("objective") = [];
coefs("acc") = [];

%%Input checks
% scenario("x_0") and scenario("x_leader"): make sure cars don't overlap
% scenario("x_0"): length == length of scenario("I_a") and scenario("I_h")
% scenario("v_0"): length == length of scenario("I_h")
% params("u_0"): length of output == length of scenario("I_a")

%% Manual optimization step

U_func = params("u_0");
U = U_func(params("t_int"));
[~, XV] = ode45(@(t,XV) F(t, XV, U_func(t)', scenario, params), params("t_int"), [scenario("x_0"); scenario("v_0")]);

V = XV(:, end - length(scenario("I_h")) + 1: end);
X = XV(:, 1:end - length(scenario("I_h")));
z = J(X, V, U_1, scenario, params)
[~, dim] = size(XV);
P_0 = zeros(dim, 1);
[t,PQ] = ode45(@(t,PQ) F_adjoint(t, PQ, X, V, U, scenario, params), flip(params("t_int")), P_0);

PQ = flip(PQ,1);
Q = PQ(:, end - length(scenario("I_h")) + 1: end);
P = PQ(:, 1:end - length(scenario("I_h")));

dz = [];

for i = scenario("I_a")
    grad = P(:, i) - L_partial(params("t_int"), i, X, V, U, "v_a", scenario, params);
    if ismember(i+1, scenario("I_h"))
        grad = grad + ACC_partial(X(:, i), X(:,i+1), U(:, scenario("I_a")==i), V(:, scenario("I_h")==i+1), 3, params) .* Q(:, scenario("I_h")==i+1); 
    end
    dz = [dz, grad];
end

figure(1)
plot(X)

U_1 = U - 0.9 * dz;
U_func = @(t) fun_interp(U_1, t);
[~, XV] = ode45(@(t,XV) F(t, XV, U_func(t)', scenario, params), params("t_int"), [scenario("x_0"); scenario("v_0")]);

V = XV(:, end - length(scenario("I_h")) + 1: end);
X = XV(:, 1:end - length(scenario("I_h")));
z = J(X, V, U_1, scenario, params)
figure(2)
plot(params("t_int"), x_leader(params("t_int")))
hold on 
plot(params("t_int"), X)
[~, dim] = size(XV);
P_0 = zeros(dim, 1);
[t,PQ] = ode45(@(t,PQ) F_adjoint(t, PQ, X, V, U_1, scenario, params), flip(params("t_int")), P_0);

PQ = flip(PQ,1);
Q = PQ(:, end - length(scenario("I_h")) + 1: end);
P = PQ(:, 1:end - length(scenario("I_h")));

%% Optimizer
options = optimoptions('fmincon','Display','iter-detailed', ...
                        'SpecifyObjectiveGradient',useGradient,...
                        'FunValCheck','on', 'DerivativeCheck', 'off',...
                        'maxfunevals',1e5,'StepTolerance',1e-10,'algorithm', 'interior-point');

fun = @(U) objective_gradient(U, scenario, params);
U_0 = params("u_0");

A = [];
b = [];
Aeq = [];
beq = [];
lb = zeros(size(U_0(params("t_int"))));
ub = 50*ones(size(U_0(params("t_int"))));
nonlcon = [];
[U_star,~,~,~,~,grad,~] = fmincon(fun, U_0(params("t_int")), A,b,Aeq,beq,lb,ub,nonlcon,options);


%% Function to optimize
function [z, dz] = objective_gradient(U, scenario, params)
    global fun_interp
    % Solve forward system 
    U_func = @(t) fun_interp(U, t);
    [~, XV] = ode45(@(t,XV) F(t, XV, U_func(t)', scenario, params), params("t_int"), [scenario("x_0"); scenario("v_0")]);
    
    V = XV(:, end - length(scenario("I_h")) + 1: end);
    X = XV(:, 1:end - length(scenario("I_h")));
    
    % Evaluate Objective 
    z = J(X, V, U_func(params("t_int")), scenario, params);
    
    % Solve backward system
    [~, dim] = size(XV);
    P_0 = zeros(dim, 1);
    [~,PQ] = ode45(@(t,PQ) F_adjoint(t, PQ, X, V, U, scenario, params), flip(params("t_int")), P_0);

    PQ = flip(PQ,1);
    Q = PQ(:, end - length(scenario("I_h")) + 1: end);
    P = PQ(:, 1:end - length(scenario("I_h")));
    
    dz = [];
    for i = scenario("I_a")
        grad = P(:, i) - L_partial(params("t_int"), i, X, V, U, "v_a", scenario, params);
        if ismember(i+1, scenario("I_h"))
            grad = grad + ACC_partial(X(:, i), X(:,i+1), U(:, scenario("I_a")==i), V(:, scenario("I_h")==i+1), 3, params) .* Q(:, scenario("I_h")==i+1); 
        end
        dz = [dz, grad];
    end 

end 


%% Forward and backward systems

function PQ_dot = F_adjoint(t, PQ, X, V, U, scenario, params) %X\in R(t_int X N+M), V \in R(t_int X M), U \in R(t_int X N)
    global fun_interp

    Q = PQ(end - length(scenario("I_h")) + 1: end);
    P = PQ(1:end - length(scenario("I_h")));
    
    P_dot = zeros(length(P), 1);
    Q_dot = zeros(length(Q), 1);
        
    % For x_h and v_h we need to know: x,v, xl, vl maybe xf,vf all as functions
    for i = scenario("I_h")
        x = fun_interp(X(:, i), t);
        v = fun_interp(V(:, scenario("I_h")==i), t);
        if i == 1 
            xl = scenario("x_leader");
            xl = xl(t);
            vl = scenario("v_leader");
            vl = vl(t); 
        elseif ismember(i - 1, scenario("I_h"))
            xl = fun_interp(X(:, i-1), t);
            vl = fun_interp(V(:, scenario("I_h")==i-1), t);
        else 
            xl = fun_interp(X(:, i-1), t);
            vl = fun_interp(U(:, scenario("I_a")==i-1), t);
        end
        P_dot(i) = L_partial(t, i, X, V, U, "x_h", scenario, params) ...
                   - ACC_partial(xl, x, vl, v, 2, params) * Q(scenario("I_h")==i);
        Q_dot(scenario("I_h")==i) = L_partial(t, i, X, V, U, "v_h", scenario, params) ...
                                    - P(i) ...
                                    - ACC_partial(xl, x, vl, v, 4, params) * Q(scenario("I_h")==i);
        if ismember(i + 1, scenario("I_h"))
            xf = fun_interp(X(:, i+1), t);
            vf = fun_interp(V(:, scenario("I_h")==i+1), t);
            P_dot(i) = P_dot(i) - ACC_partial(x, xf, v, vf, 1, params) * Q(scenario("I_h")==i+1);
            Q_dot(scenario("I_h")==i) = Q_dot(scenario("I_h")==i) - ACC_partial(x, xf, v, vf, 3, params) * Q(scenario("I_h")==i+1);        
        end 
    end 

    % For x_a need to know: maybe x,v and xf,vf all as functions
    for i = scenario("I_a") 
        P_dot(i) = L_partial(t, i, X, V, U, "x_a", scenario, params);
        if ismember(i+1, scenario("I_h"))
            x = fun_interp(X(:, i), t);
            v = fun_interp(U(:, scenario("I_a")==i), t);
            xf = fun_interp(X(:, i+1), t);
            vf = fun_interp(V(:, scenario("I_h")==i+1), t);
            P_dot(i) = P_dot(i) - ACC_partial(x, xf, v, vf, 1, params) * Q(scenario("I_h")==i+1);
        end 
    end
    
    PQ_dot = [P_dot; Q_dot];
end 

function XV_dot = F(t, XV, U, scenario, params) %X \in R(2M+N) and U \in R(N) 
    V = XV(end - length(scenario("I_h")) + 1: end);
    X = XV(1:end - length(scenario("I_h")));
    
    X_dot = zeros(length(X), 1);
    X_dot(scenario("I_h")) = V; 
    X_dot(scenario("I_a")) = U; 
    
    V_dot = zeros(length(V), 1);
    for i = scenario("I_h")
        x = X(i);
        v = V(scenario("I_h") == i);
        if i == 1
            xl = scenario("x_leader");
            xl = xl(t);
            vl = scenario("v_leader");
            vl = vl(t);
        elseif ismember(i - 1, scenario("I_h"))
            vl = V(scenario("I_h") == i-1);
            xl = X(i-1);
        else
            vl = U(scenario("I_a") == i-1);
            xl = X(i-1);
        end
        
        V_dot(scenario("I_h")==i) = ACC(xl, x, vl, v, params);
    end     
    XV_dot = [X_dot; V_dot];
end 



%% Objective function 
function j = J(X, V, U, scenario, params) % both X and U are matrices
    VU = [V, U];
    idx = scenario("I_a");
    X_diff = [];
    if idx(1) == 1
        xl = scenario("x_leader");
        X_diff = xl(params("t_int")) - X(:, 1) - params("l") - params("safe_dist");
        idx = idx(2:end);
    end 
    
    X_diff = [X_diff, X(:, idx-1) - X(:, idx) - params("l") - params("safe_dist")];
    
    integrand = (sum((VU-params("des_v")).^2, 2) + sum(1./(X_diff).^2, 2));
    
    j = diff(params("t_int"))'* integrand(1:end-1);
end 

function l_partial = L_partial(t, i, X, V, U, var, scenario, params)
    global fun_interp
    switch var
        case "x_h"
            l_partial = 0;
        
        % we need x, xl, maybe xf
        case "x_a"
            x = fun_interp(X(:, i), t);
            if i == 1
                xl = scenario("x_leader");
                xl = xl(t);
            else 
                xl = fun_interp(X(:, i-1), t);
            end 
            l_partial = 2./(xl - x - params("l") - params("safe_dist")).^3;
            if ismember(i+1, scenario("I_a"))
                xf = fun_interp(X(:, i+1), t);
                l_partial = l_partial - 2./(x - xf - params("l") - params("safe_dist")).^3;
            end 
 
        case "v_h"
            v = fun_interp(V(:, scenario("I_h")==i), t);
            l_partial = 2*(v - params("des_v"));
            
        case "v_a"
            v = fun_interp(U(:, scenario("I_a")==i), t);
            l_partial = 2*(v - params("des_v"));
    end 
end 


%% Car following model 
function a = ACC(xl, xf, vl, vf, params)
    V = @(head_way) params("v_max") * ((tanh(head_way - params("safe_dist")) + tanh(params("l") + params("safe_dist")))./(1 + tanh(params("l") + params("safe_dist"))));
    a = params("alpha") *(V(xl - xf - params("l")) - vf) + params("beta") * (vl - vf) ./ ((xl - xf - params("l")).^2);
end 

function a_partial = ACC_partial(xl, xf, vl, vf, var, params)
    
    V_prime = @(head_way) (params("v_max") * (sech(head_way - params("safe_dist"))).^2) ./ (1+tanh(params("l")+params("safe_dist")));
    head_way = xl - xf - params("l");
    switch var
        case 1 
            a_partial = params("alpha") * V_prime(head_way) - 2*params("beta")*((vl-vf)./head_way.^3);
        case 2
            a_partial = -1*params("alpha") * V_prime(head_way) + 2*params("beta")*((vl-vf)./head_way.^3);
        case 3
            a_partial = params("beta") ./ head_way.^2;
        case 4
            a_partial = -1*params("alpha") - params("beta") ./ head_way.^2;
    end 
end 


