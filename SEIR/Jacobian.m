function J = Jacobian(para,y1,y3)
% Jacobian matrix \partial f/\partial x
J = zeros(4);
beta = para(1);
sigma = para(2);
gamma = para(3);

 S = y1;
 I = y3;
J(1,1) = -beta * I;
J(1,3) = -beta * S;

J(2,1) = beta * I;
J(2,2) = -sigma;
J(2,3) = beta * S;

J(3,2) = sigma;
J(3,3) = -gamma;

J(4,3) = gamma;
end