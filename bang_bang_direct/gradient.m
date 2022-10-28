function grad = gradient(z, auxdata)
%OBJECTIVE Numerical solution with direct method, evaluates velocities
%between time frames

    %Extract variables from z setup
    N = auxdata.N;
    T = auxdata.T_size;
    nvars = auxdata.nvars
    
    start_v = N+1;
    start_u = start_v + N+1;

    x = z(1:N+1);
    v = z(start_v + (1:N+1));
    u = z(start_u + (1:N));

    % TO-DO: Calculate gradient
    grad                = zeros(nvars,1) ;
    grad(start_v+1)     = -1/2 ;
    grad(start_v+(2:N)) = -1   ;
    grad(start_v+N+1)   = -1/2 ;
    
end