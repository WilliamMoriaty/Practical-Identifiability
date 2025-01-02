function S = Sens_NN(alpha,beta,omega,N)
% Jacobian matrix of NN
% data size: N
% x\in [0,1]

M = size(alpha,1);
S_alpha = zeros(N,M);
S_beta = zeros(N,M);
S_omega = zeros(N,M);

for i = 1:N
    for j= 1:M
        S_alpha (i,j) = omega(j) * (1-(tanh(alpha(j)*((i-1)/(N-1))+beta(j)))^2) * (i-1)/(N-1);
        S_beta (i,j) = omega(j) * (1-(tanh(alpha(j)*((i-1)/(N-1))+beta(j)))^2);
        S_omega (i,j) = tanh(alpha(j)*((i-1)/(N-1))+beta(j));        
    end
end

S = [S_alpha S_beta S_omega];
end


