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
        S_alpha (i,j) = omega(j) * Dsigma(alpha(j)*((i-1)/(N-1))+beta(j)) * (i-1)/(N-1);
        S_beta (i,j) = omega(j) * Dsigma(alpha(j)*((i-1)/(N-1))+beta(j));
        S_omega (i,j) = sigma(alpha(j)*((i-1)/(N-1))+beta(j));        
    end
end

S = [S_alpha S_beta S_omega];
end

function Ds = Dsigma(x)
% Derivative of Relu(x)
if x > 0
    Ds = 1;
else
    Ds = 0;
end
end

function s = sigma(x)
if x > 0
    s = x;
else
    s = 0;
end
end