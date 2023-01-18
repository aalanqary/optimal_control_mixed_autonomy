function f = objective(U, auxdata)  
    [time_v, v] = system_solve(U, auxdata);
    
    figure(1)
    plot(time_v, v)
    drawnow
    figure(2)
     plot(auxdata.tau, U)
     plot(U)
     drawnow
     pause(3)
    f = -trapz(time_v, v);
    v = griddedInterpolant(time_v, v, "previous");
    v = v(auxdata.tau);
    const(U, auxdata)
    %assert(sum(abs(U) <= auxdata.g + auxdata.k3 * v.^2) == length(v)); 
    if sum(abs(U) <= auxdata.g + auxdata.k3 * v.^2) ~= length(v)
        disp(abs(U) <= auxdata.g + auxdata.k3 * v.^2);
    end
end

