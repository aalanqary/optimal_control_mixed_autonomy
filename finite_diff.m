function grad = finite_diff(func, X,  h)
    grad = 0 * X; 
    for j = 1:1:2
        for i = 1:1:length(X)
            X_plus = X; 
            X_minus = X; 
            X_plus(i, j) = X(i, j) + h; 
            X_minus(i, j) = X(i, j) - h; 
            f_plus = func(X_plus); 
            f_minus = func(X_minus); 
            grad(i, j) = (f_plus - f_minus) / (2*h); 
        end 
    end 
end 
