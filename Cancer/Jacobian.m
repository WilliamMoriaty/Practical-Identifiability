function J = Jacobian(para,E,T,C)
% Jacobian matrix \partial f/\partial x
J = zeros(3);
sigma = para(1);rho = para(2);eta = para(3);
mu = para(4);eps = para(5);
beta1 = para(6);beta2 = para(7);
phi = para(8);lambda = para(9);
psi = para(10);
J(1,1) = -sigma - mu * T;
J(1,2) = -rho * C./((eta+T)^2) - mu*E;
J(1,3) = rho./(eta+T) + eps;

J(2,1) = - phi * T;
J(2,2) = (beta1-2*beta1*beta2*T) - phi*E;
J(2,3) = lambda;

J(3,1) = mu * T;
J(3,2) = mu * E;
J(3,3) = -psi;

end