v00 = @(t) 0*t+0.25;
useGradient = true;
t0 = 0;
t1 =10;
nt = 70;
t_int = linspace(t0, t1, nt); 

%% Define the optimization problem 
options = optimoptions('fmincon','Display','iter-detailed', ...
                        'SpecifyObjectiveGradient',useGradient,...
                        'FunValCheck','on', ...
                        'maxfunevals',1e5,'StepTolerance',1e-16);
if useGradient
    fun = @(v) objective_gradient(v,t_int,t1);
else
    fun = @(v) objective(v,t_int,t1);
end 

A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon = [];
[v2,~,~,~,~,grad,~] = fmincon(fun,v00(t_int),A,b,Aeq,beq,0*ones(1,nt),50*ones(1,nt),nonlcon, options);

%% Solve the system with the optimal velocity 
x0 = [5; 0];
v0 = [0; 0];

vl =@(t) interp1(t_int, v2, t);

[t2,y]=ode45(@(t,y) vl(t),t_int,x0(1));
xl=@(t) interp1(t2, y(:,1), t);

alpha=1;
beta=1;
l=4;
v_max = 14;
safe_des = 3;

V=@(x)v_max * ((tanh(x- safe_des) + tanh(x + safe_des))./(1 + tanh(l + safe_des)));
ACC=@(c,d,t) alpha*(V(xl(t)-c-l)-d)+beta*(vl(t)-d)./((xl(t)-c-l).^2);
[t,y] = ode45(@(t,y)[y(2);ACC(y(1),y(2),t)],t2,[x0(2);v0(2)]);


%% Plot results 
figure(1)
plot(t,xl(t),t,y(:,1))
figure(2)
plot(t,vl(t),t,y(:,2))
figure(3)
plot(t,vl(t))

%% Define functions 
function z = objective(v,t_int,t1)
x0 = [5; 0];
v0 = [0; 0];

vl =@(t) interp1(t_int, v, t);


[t2,y]=ode45(@(t,y) vl(t), t_int,x0(1));

xl=@(t) interp1(t2, y(:,1), t);


alpha=1;
beta=1;
l=4;
v_max = 14;
safe_des = 3;

V=@(x) v_max * ((tanh(x - safe_des) + tanh(x + safe_des))./(1 + tanh(l + safe_des)));
ACC=@(c,d,t) alpha*(V(xl(t)-c-l)-d)+beta*(vl(t)-d)./((xl(t)-c-l).^2);

[t,y]=ode45(@(t,y)[y(2);ACC(y(1),y(2),t)],t2,[x0(2);v0(2)]);


des_vel = 3;
z = t1*sum(diff(t).*(y(1:end-1,2)-des_vel).^2);%+ %0.005*t1*sum(diff(t).*(ACC(y(1:end-1,1),y(1:end-1,2),t(1:end-1))).^2)
end


function [z, dz] = objective_gradient(v,t_int,t1)
x0 = [5; 0];
v0 = [0; 0];

vl =@(t) interp1(t_int, v, t);


[t2,y]=ode45(@(t,y) vl(t), t_int,x0(1));

xl = @(t) interp1(t2, y(:,1), t);


alpha=1;
beta=1;
l=4;
v_max = 14;
safe_des = 3;

V = @(x) v_max * ((tanh(x - safe_des) + tanh(x + safe_des))./(1 + tanh(l + safe_des)));
ACC=@(c,d,t) alpha*(V(xl(t)-c-l)-d)+beta*(vl(t)-d)./((xl(t)-c-l).^2);

[t3,y] = ode45(@(t,y)[y(2);ACC(y(1),y(2),t)],t2,[x0(2);v0(2)]);

des_vel = 3;
z = t1*sum(diff(t3).*(y(1:end-1,2)-des_vel).^2); %+ %0.005*t1*sum(diff(t).*(ACC(y(1:end-1,1),y(1:end-1,2),t(1:end-1))).^2)


xf = @(t) interp1(t3, y(:,1), t);
vf = @(t) interp1(t3, y(:,2), t);

V_prime = @(x) (v_max * (sech(x - safe_des)).^2) / (1+tanh(l+safe_des));

f_a = @(t,p) [-1*(alpha*V_prime(xl(t)-xf(t) - l) - beta*(2*(vl(t) - vf(t))/(xl(t)-xf(t) - l).^3))*p(3);
    (alpha*V_prime(xl(t)-xf(t) - l) - beta*(2*(vl(t) - vf(t))/(xl(t)-xf(t) - l).^3))*p(3);
     2*(vf(t) - des_vel) - p(2) + (alpha + beta/(xl(t)- xf(t) - l).^2)*p(3)];

p0 = [0; 0; 0];

[t4,p] = ode45(f_a, t2, p0);

dz =  p(:, 1) + ((beta ./ (xl(t4) - xf(t4) - l).^2) .* p(:, 3));

end

