function y = hill_function(x, Vmax, Kd, n)
    

    % 
    y = (Vmax * x.^n) ./ (Kd^n + x.^n);
end
