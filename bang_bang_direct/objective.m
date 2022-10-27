function f = objective(z)
%OBJECTIVE Numerical solution with direct method, evaluates velocities
%between time frames

    %Extract variables from z setup
    N = auxdata.N;
    T = auxdata.T_size;
    
    start_v = N+1;
    start_u = start_v + N+1;

    x = z(1:N+1);
    v = z(start_v + (1:N+1));
    u = z(start_u + (1:N));

    % Velocities calculated between intervals due to convention
    vave = (v(2:end) + v(1:end-1))/2;

    f = -sum(vave);
end

