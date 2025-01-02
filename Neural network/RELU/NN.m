function f = NN(x,alpha,beta,omega)
% Neural Network
M = size(alpha,1);
f = 0;
for j = 1:M
f = f + omega(j) * sigma(alpha(j) * x + beta(j));
end

end

function s = sigma(x)
if x > 0
    s = x;
else
    s = 0;
end
end