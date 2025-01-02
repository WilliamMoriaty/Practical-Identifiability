function f = NN(x,alpha,beta,omega)
% Neural Network
M = size(alpha,1);
f = 0;
for j = 1:M
f = f + omega(j) * tanh(alpha(j) * x + beta(j));
end

end

