function f = plot_x_v_u(auxdata, z)
    %Extract variables from z setup
    N = auxdata.N;
    T = auxdata.T_size;
    
    start_v = N+1;
    start_u = start_v + N+1;

    x = z(1:N+1);
    v = z(start_v + (1:N+1));
    u = z(start_u + (1:N));

    axis = linspace(0, N+1)

    disp(N)

    plot(x);
    figure;
    plot(v)
    figure;
    plot(u)