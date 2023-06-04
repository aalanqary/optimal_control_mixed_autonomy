%% Define problem 

% Platoon params
    auxdata.platoon = [1, 1, 1, 1];
    auxdata.len_platoon = length(auxdata.platoon);
    auxdata.Ia = find(auxdata.platoon);
    auxdata.Ih = find(auxdata.platoon - 1);

% Bando-FtL params
    auxdata.safe_dist = 4; 
    auxdata.v_max = 35; 
    auxdata.alpha = 0.1;
    auxdata.beta = 21*25;
    auxdata.k = 0.2;
    auxdata.l = 5;

% Optimized Bando-FtL params
%     auxdata.safe_dist = 3.92988471e+01; 
%     auxdata.v_max = 40; 
%     auxdata.alpha = 5.64872416e-04;
%     auxdata.beta = 2.33052191e+02;
%     auxdata.k = 1.98761985e-01;
%     auxdata.l = 5; 

% objective function params
    auxdata.mu1 = 1000;
    auxdata.mu2 = 0;
    auxdata.mu3 = 0;

% Energy function params
%     load(['RAV4_coeffs.mat']);
%     auxdata.C = [C0, C1, double(C2), C3];
%     auxdata.p = [p0, p1, p2];
%     auxdata.q = [double(q0), q1];

%Arctan Barrier auxdata (a(-arctan(bx+c)+pi/2)
    auxdata.a = 1;
    auxdata.b = 50;
    auxdata.c = 0.6;
    auxdata.min_translation = 3.5;
    auxdata.c = -auxdata.b*auxdata.min_translation + auxdata.c;

%Arctan max velocity auxdata (a(-arctan(bx+c)+pi/2))
    auxdata.d = 1;
    auxdata.e = 3;
    auxdata.f = -12.5;
    auxdata.max_translation = 300;
    auxdata.f = -auxdata.e*auxdata.max_translation + auxdata.f;

% const grad params
    auxdata.eps = 0.1;
    auxdata.gamma = 0.1;
    auxdata.d_min = 2.5;
    auxdata.d_max = 300;

% Basic Leader traj
%     
%     T = 100;
%     auxdata.dt = 0.1;
%     auxdata.udt = 1; 
%     auxdata.utime = (0:auxdata.udt:T)';
%     auxdata.time = (0:auxdata.dt:T)';
%     vl = @(t) (t<=1000).*20;
%     vl = vl(auxdata.time); 
%     vl = smoothdata(vl, "movmean", 50);
%     auxdata.vl = griddedInterpolant(auxdata.time, vl);

% %Simple leader traj
     T = 250; 
     auxdata.dt = 0.1;
     auxdata.udt = 1; 
     auxdata.utime = (0:auxdata.udt:T)';
     auxdata.time = (0:auxdata.dt:T)';
     
     vl = @(t) (t<=80) .* 30 ...
                 + (((t>80) + (t<= 120)) ==2) .* (-t./8 + 40) ...
                 + (((t>120) + (t<= 150)) ==2) .* 25 ...
                 + (((t>150) + (t<= 190)) ==2) .* (t./8 + 25/4) ...
                 + (t>190) .* 30;
     vl = vl(auxdata.time); 
     vl = smoothdata(vl, "movmean", 50);
     auxdata.vl = griddedInterpolant(auxdata.time, vl);

%Sinusoidal leader traj
%     T = 250; 
%     auxdata.dt = 0.1;
%     auxdata.udt = 1; 
%     auxdata.utime = (0:auxdata.udt:T)';
%     auxdata.time = (0:auxdata.dt:T)';
%     
%     auxdata.vl = @(t) sin(0.1*t) + cos(0.05*t)  - t/100 + 32;

%Real leader traj
%     T = 704; 
%     auxdata.dt = 0.1;
%     auxdata.udt = 1; 
%     auxdata.utime = (0:auxdata.udt:T)';
%     auxdata.time = (0:auxdata.dt:T)';
%     
%     data = readtable("data_v2_preprocessed_west/2021-04-22-12-47-13_2T3MWRFVXLW056972_masterArray_0_7050.csv");
%     vl = data.Velocity * (1000/3600); 
%     vl = smoothdata(vl,'movmean',200);
%     auxdata.vl = griddedInterpolant(auxdata.time,vl(1:length(auxdata.time)));

% Initial conditions and leader position 
    eq = round(eq_headway(auxdata.vl(0), auxdata), 5);
    x0 = (eq+auxdata.l) * flip(0:1:auxdata.len_platoon)';
    % x0 = 1*auxdata.vl(0) * flip(0:1:auxdata.len_platoon)';
    x0 = x0 - x0(2); 
    auxdata.v0 = ones(auxdata.len_platoon, 1) * auxdata.vl(0); 
    opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
    [~, xl] = ode45(@(t,x) auxdata.vl(t), auxdata.time, 0, opts);
    xl = xl + x0(1);
    auxdata.x0 = x0(2:end);
    auxdata.xl = griddedInterpolant(auxdata.time, xl);

%Initial Guess U0  
    U0 = diff(auxdata.vl(auxdata.utime)) ./ auxdata.udt;
    U0 = [U0;0];
    U0 = repmat(U0, [1, length(auxdata.Ia)]);
    [X0, V0, A0] = generic_system_solve(U0, auxdata);

%% Run Optimizaer 

options = optimoptions('fmincon','Display','iter-detailed', ...
                        'SpecifyObjectiveGradient', true ,...
                        'SpecifyConstraintGradient',true, ...
                        'GradConstr', 'on',...
                        'CheckGradients',false,...
                        'FunValCheck','off', 'DerivativeCheck', 'off',...
                        'maxfunevals',1e6, 'StepTolerance',1e-8, ...
                        'algorithm', 'sqp', ...
                        'ConstraintTolerance', 1e-10);

nonlcon = @(U) nonlcon_const_max(U, auxdata);
fun = @(U) objective_gradient_acc(U, auxdata); 
% nonlcon = []; %@(U) nlc(U, auxdata);
A = [];
b = [];
% [A, b] = lc_v(auxd
% ata); 
Aeq = []; beq = []; 
a_min = []; %-3 * ones(size(U0)); 
a_max = []; %3 * ones(size(U0));
tic
[U_star, f_val, ~, output, ~, grad] = fmincon(fun, U0, A, b, Aeq, beq, a_min, a_max, nonlcon, options); % switch to ip solver instead of sqp
[X_star, V_star, A_star] = generic_system_solve(U_star, auxdata);
toc

