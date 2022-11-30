function f = objective(U, auxdata)  
    [time_v, v] = system_solve(U, auxdata);
%     figure(1)
%     plot(time_v, v)
%     drawnow
%     figure(2)
%     plot(auxdata.tau, U)
%     drawnow
    f = -trapz(time_v, v);
end

