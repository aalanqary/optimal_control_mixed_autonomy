function y = Euler(h, y_0, t_int)
t_int = t_int(1):h:t_int(2);
y = zeros(size(t_int));  % allocate the result y
y(1) = y_0;  % the initial y value
n = numel(y);  % the number of y values
% The loop to solve the DE
for i=1:n-1
    f = -15 * y(i); 
    y(i+1) = y(i) + h * f;
end