function [dc, dceq] = const_eq_grad(U, auxdata)  
    [time_v, v] = system_solve(U, auxdata);
    v = griddedInterpolant(time_v, v, "linear");
    U = griddedInterpolant(auxdata.tau(1:end-1), U(1:end-1), "linear");
    df = zeros(3, auxdata.N); 

    %For first constraint: 
    lambdaT_1 = 1;
    Lambda_1 = ode5(@(t,lambda) lambda *(auxdata.k1 + 2*auxdata.k2*v(t)), flip(time_v), lambdaT_1); 
    Lambda_1 = griddedInterpolant(time_v, Lambda_1, "linear");
    df1 = arrayfun(@(a,b) integral(@(t) Lambda_1(t), a, b), auxdata.tau(1:end-1), auxdata.tau(2:end));  
    if isrow(df1)
        df1 = df1';
    end 
    dceq = [df1; 0];
    
    dc = [];
end