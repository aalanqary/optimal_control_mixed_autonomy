function [X, V, A] = system_solve(U_vec, auxdata) 
%%IMPORTANT NOTE: this only solve systems that have AV in index 1 then all HV
    
    if isempty(U_vec)
        U = @(t) [];
    else
        U = griddedInterpolant(auxdata.utime, U_vec, "previous");
    end

    
    F = @(t, XV) [XV(auxdata.len_platoon+1:end); ...
        U(t); ACC(XV(auxdata.Ih - 1), XV(auxdata.Ih), XV(auxdata.len_platoon + auxdata.Ih - 1), XV(auxdata.len_platoon + auxdata.Ih), auxdata)];
   

    XV = ode5(@(t, XV) [XV(auxdata.len_platoon+1:end); U(t); ACC(XV(auxdata.Ih - 1), XV(auxdata.Ih), XV(auxdata.len_platoon + auxdata.Ih - 1), XV(auxdata.len_platoon + auxdata.Ih), auxdata)], ...
                        auxdata.time, [auxdata.x0; auxdata.v0]);
    V = XV(:, end - auxdata.len_platoon + 1: end);
    X = XV(:, 1:end - auxdata.len_platoon);

    if nargout > 2
        A = zeros(size(X));
        Xl = [auxdata.xl(auxdata.time), X]; 
        Vl = [auxdata.vl(auxdata.time), V];  
        Ah = ACC(Xl(:, auxdata.Ih), X(:, auxdata.Ih), ... 
                 Vl(:, auxdata.Ih), V(:, auxdata.Ih), auxdata); 
        A(:, auxdata.Ih) = Ah;
        if ~isempty(U)
            A(:, auxdata.Ia) = U(auxdata.time); 
        end 
    end 

