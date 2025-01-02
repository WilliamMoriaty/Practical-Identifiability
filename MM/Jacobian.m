function J = Jacobian(para,y1,y2)
% Jacobian matrix \partial f/\partial x
J = zeros(4);
k1 = para(1);
k2 = para(2);
k3 = para(3);
J(1,1) = -k1 * y2;
J(1,2) = -k1 * y1;
J(1,3) = k2;

J(2,1) = -k1 * y2;
J(2,2) = -k1 * y1;
J(2,3) = k2 + k3;

J(3,1) = k1 * y2;
J(3,2) = k1 * y1;
J(3,3) = -(k2 + k3);

J(4,3) = k3;
end