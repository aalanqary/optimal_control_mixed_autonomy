function df = obj_grad(U, auxdata)  
    [time_v, v] = system_solve(U, auxdata);
    lambdaT = 0;
    v = griddedInterpolant(time_v, v, "linear");
    Lambda = ode5(@(t,lambda) 1 + lambda *(auxdata.k1 + 2*auxdata.k2*v(t)), flip(time_v), lambdaT);    
    Lambda = flip(Lambda);
    Lambda = griddedInterpolant(time_v, Lambda, "linear");
    df = arrayfun(@(a,b) integral(@(t) Lambda(t), a, b), auxdata.tau(1:end-1), auxdata.tau(2:end));  
    df = [df, 0];
end
