function dxdthe2 = partial_theta2(t,x,para)
    % \partial x/\partial \theta_1
    %[~, Y] = ode45(@(t,y)SEIR_eq(t,y,para), tspan1, y0);
    %y1_t = Y(end,1);y2_t = Y(end,2);y3_t = Y(end,3); y4_t = Y(end,4);
    %J = Jacobian(para,y1_t,y3_t);
     beta = para(1);
sigma = para(2);
gamma = para(3);

    S = x(5);
    E = x(6);
    I = x(7);
    % ODE
    dxdthe2 = zeros(8,1);
    J = Jacobian(para,S,I);
    dxdthe2(1) = J(1,1)*x(1) + J(1,2)*x(2) + J(1,3)*x(3) + J(1,4)*x(4);
    dxdthe2(2) = -E + J(2,1)*x(1) + J(2,2)*x(2) + J(2,3)*x(3) + J(2,4)*x(4);
    dxdthe2(3) = E + J(3,1)*x(1) + J(3,2)*x(2) + J(3,3)*x(3) + J(3,4)*x(4);
    dxdthe2(4) = 0 + J(4,1)*x(1) + J(4,2)*x(2) + J(4,3)*x(3) + J(4,4)*x(4);
    dxdthe2(5) = -beta * S * I ;  % dy1/dt
    dxdthe2(6) = beta * S * I  - sigma * E; % dy2/dt
    dxdthe2(7) = sigma * E - gamma * I; % dy3/dt
    dxdthe2(8) = gamma * I; % dy4/dt
end