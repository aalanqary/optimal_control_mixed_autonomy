global T N x0 x1

T = 1; N = 10; x0 = [0; 0]; x1 = [0; pi/3];

u_init = zeros(1, N);

u_opt = fmincon(@cost3, u_init, [], [], [], [], [], [], @nonlcon);
% u_opt = griddedInterpolant(linspace(0, T, N), u_opt, "next");
[t, x] = ode45(@my_ode2, [0, T], x0, [], u_opt);

function xdot = my_ode2(t, x, u)
    uu = udet2(t, u);
%     uu = u(t);
    xdot = zeros(2, 1);
    xdot(1) = uu;
    xdot(2) = sqrt(1+uu^2);
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
%     u = griddedInterpolant(linspace(0, T, N), u, "next");
    [t, x] = ode45(@my_ode2, [0, T], x0, [], u);
    J2 = trapz(t, x(:, 1));
end

function [c, ceq] = nonlcon(u)
    global T x0 N
%     u = griddedInterpolant(linspace(0, T, N), u, "next");
    [t, x] = ode45(@my_ode2, [0, T], x0, [], u);
%     n = length(x);
    x1 = x(:, 1)';
    x2 = x(:, 2)';
    c = [];
    ceq(1) = x1(end);
    ceq(2) = x2(end) - pi/3;
    ceq = ceq';
end 