function J = Jacobian(para,X,Y)
% Jacobian matrix \partial f/\partial x
J = zeros(2);
J(1,1) = para(1) - para(2) * Y;
J(1,2) = -para(2) * X;
J(2,1) = para(3) * Y;
J(2,2) = para(3) * X - para(4);

end