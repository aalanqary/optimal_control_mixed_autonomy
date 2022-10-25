function f = eq_const(U, auxdata)

    [X, V] = system_solve(U, params);

  % extract aux parameters
  g  = auxdata.g ;
  h  = auxdata.h ;
  k0 = auxdata.k0 ;
  k1 = auxdata.k1 ;
  k2 = auxdata.k2 ;
  k3 = auxdata.k3 ;

  vave   = (v(2:end)+v(1:end-1))/2 ;
  deltax = x(2:end)-x(1:end-1) ;
  deltav = v(2:end)-v(1:end-1) ;

  f = [ deltax - h*vave ; ...
        deltav - h*(u - k0 - vave.*(k1+k2*vave)) ; ...
        x(1) ; v(1) ; v(end) ; ...
        (g+k3*vave.^2)+u ; ...
        (g+k3*vave.^2)-u ] ;
end