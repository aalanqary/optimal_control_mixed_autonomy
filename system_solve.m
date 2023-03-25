function [X, V, A] = system_solve(U_vec, auxdata) 
%%IMPORTANT NOTE: this only solve systems that have AV in index 1 then all
%%HV and time spacing of AV is 1 s
    
    if isempty(U_vec)
        U = @(t) []* t;
    else
        U = griddedInterpolant(auxdata.utime, U_vec, "previous");
    end
    Ih = auxdata.Ih;
    len_platoon = auxdata.len_platoon;
    time = auxdata.time;
    if ~isempty(U_vec)
        F = @(t, XV) [XV(len_platoon+1:end); ...
            U_vec(floor(t)+1); ACC(XV(Ih - 1), XV(Ih), XV(len_platoon + Ih - 1), XV(len_platoon + Ih), auxdata)];
    else 
        F = @(t, XV) [XV(len_platoon+1:end); ...
            ACC(XV(Ih - 1), XV(Ih), XV(len_platoon + Ih - 1), XV(len_platoon + Ih), auxdata)];
    end 

    XV = ode3(F, auxdata.time, [auxdata.x0; auxdata.v0]);
    V = XV(:, end - len_platoon + 1: end);
    X = XV(:, 1:end - len_platoon);

    if nargout > 2
        A = zeros(size(X));
        Xl = [auxdata.xl(time), X]; 
        Vl = [auxdata.vl(time), V];  
        Ah = ACC(Xl(:, Ih), X(:, Ih), ... 
                 Vl(:, Ih), V(:, Ih), auxdata); 
        A(:, Ih) = Ah;
        if ~isempty(U_vec)
%             U_tall = repmat(U_vec, 1, 10)';
            A(:, auxdata.Ia) = U(time); 
        end 
    end 

