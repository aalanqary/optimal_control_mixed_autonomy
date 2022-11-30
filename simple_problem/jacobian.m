function jac = jacobian(U, auxdata)  
    [time_v, v] = system_solve(U, auxdata);
    v = griddedInterpolant(time_v, v, "linear");
    U = griddedInterpolant(auxdata.tau(1:end-1), U(1:end-1), "linear");

    %For first constraint: 
    lambdaT_1 = 1;
    Lambda_1 = ode5(@(t,lambda) lambda *(auxdata.k1 + 2*auxdata.k2*v(t)), flip(time_v), lambdaT_1); 
    Lambda_1 = griddedInterpolant(time_v, flip(Lambda_1), "linear");
    df1 = arrayfun(@(a,b) integral(@(t) Lambda_1(t), a, b), auxdata.tau(1:end-1), auxdata.tau(2:end));  
    if iscolumn(df1)
        df1 = df1';
    end 
    df1 = [df1, 0];
    
    %For second constraint: 
    lambdaT_2 = 0;
    Lambda_2 = ode5(@(t,lambda) f_adj_2(U(t), v(t), lambda, auxdata), flip(time_v), lambdaT_2); 
    Lambda_2 = griddedInterpolant(time_v, flip(Lambda_2), "linear");
    df2 = arrayfun(@(a,b) integral(@(t) -f_int_2(U(t), v(t), Lambda_2(t), auxdata), a, b), auxdata.tau(1:end-1), auxdata.tau(2:end));  
    if iscolumn(df2)
        df2 = df2';
    end 
    df2 = [df2, 0];

    %For third constraint: 
    lambdaT_3 = 0;
    Lambda_3 = ode5(@(t,lambda) f_adj_3(U(t), v(t), lambda, auxdata), flip(time_v), lambdaT_3); 
    Lambda_3 = griddedInterpolant(time_v, flip(Lambda_3), "linear");
    df3 = arrayfun(@(a,b) integral(@(t) -f_int_3(U(t), v(t), Lambda_3(t), auxdata), a, b), auxdata.tau(1:end-1), auxdata.tau(2:end));  
    if iscolumn(df3)
        df3 = df3';
    end 
    df3 = [df3, 0];
    
    % should I just concatenate like this
    jac = sparse([df1; df2; df3]);
end

function f = f_adj_2(u, v, lambda, auxdata)
    h = auxdata.g + auxdata.k3*v^2 - u; 
    if h < -auxdata.eps
        f = 2*auxdata.k3*v;
    elseif (h >= -auxdata.eps) && (h<= auxdata.eps)
        f = (-1/4*auxdata.eps) * (4*auxdata.k3^2 * v^3 + 4*auxdata.k3*v*(auxdata.g - u - auxdata.eps));
    else
        f = 0; 
    end
    f = f + lambda * (auxdata.k1 + 2*auxdata.k2 * v);
end 

function f = f_adj_3(u, v, lambda, auxdata)
    h = auxdata.g + auxdata.k3*v^2 + u; 
    if h < -auxdata.eps
        f = 2*auxdata.k3*v;
    elseif (h >= -auxdata.eps) && (h<= auxdata.eps)
        f = (-1/4*auxdata.eps) * (4*auxdata.k3^2 * v^3 + 4*auxdata.k3*v*(auxdata.g + u - auxdata.eps));
    else
        f = 0; 
    end
    f = f + lambda * (auxdata.k1 + 2*auxdata.k2 * v);
end 

function f = f_int_2(u, v, lambda, auxdata)
    h = auxdata.g + auxdata.k3*v.^2 - u; 
    bound_h = (-auxdata.eps <= h) + (auxdata.eps >= h);
    
    f = -1 * (h < -auxdata.eps) + ...
        (-(1/4*auxdata.eps) * (2 * u - 2*(auxdata.g + auxdata.k3*v.^2 - auxdata.eps))) .* (bound_h >= 2); 
    f = -f + lambda; 
end 


function f = f_int_3(u, v, lambda, auxdata)
    h = auxdata.g + auxdata.k3*v.^2 + u; 
    bound_h = (-auxdata.eps <= h) + (auxdata.eps >= h);
    f = 1 * (h < -auxdata.eps) + ...
        (-(1/4*auxdata.eps) * (2 * u + 2*(auxdata.g + auxdata.k3*v.^2 - auxdata.eps))) .* (bound_h >= 2); 
    f = -f + lambda;
end 