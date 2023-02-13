function f = objective(U, auxdata)  
    [time_v, v] = system_solve(U, auxdata);
    f = -trapz(time_v, v);
    
end

