function [X, V] = ode1_sys_solve(U_vec, auxdata)
    if isempty(U_vec)
        U = @(t) []* t;
    else
        U = griddedInterpolant(auxdata.utime, U_vec, "previous");
    end
    
    if ~isempty(U_vec)
        F = @(t, XV) [XV(auxdata.len_platoon+1:end); ...
            U(t); ACC(XV(auxdata.Ih - 1), XV(auxdata.Ih), XV(auxdata.len_platoon + auxdata.Ih - 1), XV(auxdata.len_platoon + auxdata.Ih), auxdata)];
    else 
        F = @(t, XV) [XV(auxdata.len_platoon+1:end); ...
            ACC(XV(auxdata.Ih - 1), XV(auxdata.Ih), XV(auxdata.len_platoon + auxdata.Ih - 1), XV(auxdata.len_platoon + auxdata.Ih), auxdata)];
    end 

    XV = ode3(F, auxdata.utime, [auxdata.x0; auxdata.v0]);
    V = XV(:, end - auxdata.len_platoon + 1: end);
    X = XV(:, 1:end - auxdata.len_platoon);
end 