function f = constraints(z, auxdata)
   
    %Extract variables from z setup
    N = auxdata.N;
    T = auxdata.T_size;
    
    start_v = N+1;
    start_u = start_v + N+1;

    x = z(1:N+1);
    v = z(start_v + (1:N+1));
    u = z(start_u + (1:N));

    % Define variables g through k3
    g  = auxdata.g ;
    h  = auxdata.h ;
    k0 = auxdata.k0 ;
    k1 = auxdata.k1 ;
    k2 = auxdata.k2 ;
    k3 = auxdata.k3 ;

    vave = (v(2:end)+v(1:end-1))/2 ;
    delta_x = x(2:end) - x(1:end-1);
    delta_v = v(2:end) - v(1:end-1);

    % Constraints extracted from page 159 of Numerical Optimal Controls
    % paper. Discretized as follow
    % Using uk+1/2 potentially because inputs are
    c1 = delta_x + h*vave;
    c2 = delta_v - h*(u - k0 - vave.*(k1+k2*vave)) ; 
    c3 = x(1);
    c4 = v(1);
    c5 = v(end);
    c6 = (g + k3*vave.^2) + u;
    c7 = (g+k3*vave.^2)-u;

    % Return constraints
    f = [c1; c2; c3; c4; c5; c6; c7] ;
   
end