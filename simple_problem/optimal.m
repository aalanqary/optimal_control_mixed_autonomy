% 
function opt = optimal(high, low, time, fin)
    opt = ones(1,time)*high;
    opt(time+1:fin) = low;
end