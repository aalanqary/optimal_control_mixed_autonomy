global T N x0 x1 k0 k1 k2

T = 1; N = 10; x0 = [0; 0]; x1 = [0; pi/3];
g = 9.81;
k0 = g*0.02;
k1 = g*1e-5;
k2 = g*1e-4;

u_init = zeros(1, N);

u_opt = fmincon(@cost3, u_init, [], [], [], [], -4*ones(1, N), 4*ones(1, N), @nonlcon);
u_opt = griddedInterpolant(linspace(0, T, N), u_opt, "previous");
[t, x] = ode45(@my_ode2, [0, T], x0, [], u_opt);

function xdot = my_ode2(t, x, u)
    global k0 k1 k2
%     uu = udet2(t, u);
    uu = u(t);
    xdot = zeros(2, 1);
    xdot(1) = 0;
    xdot(2) = uu - k0 - k1*x(2) - k2*x(2)^2;
end 

function vectu = udet2(t, u)
    global T N 
    a = N*t/T;
    if a == 0
        vectu = u(:, 1);
    else
        vectu = u(:, abs(ceil(a)));
    end 
end

function J2=cost3(u)
    global T x0 N
    u = griddedInterpolant(linspace(0, T, N), u, "previous");
    [t, x] = ode45(@my_ode2, [0, T], x0, [], u);
    J2 = trapz(t, -x(:, 2));
end

function [c, ceq] = nonlcon(u)
    global T x0 N
    u = griddedInterpolant(linspace(0, T, N), u, "previous");
    [t, x] = ode45(@my_ode2, [0, T], x0, [], u);
%     n = length(x);
    x1 = x(:, 1)';
    x2 = x(:, 2)';
    c = [];
    ceq(1) = x2(end);
%     ceq(2) = x2(end) - pi/3;
%     ceq = ceq';
end 