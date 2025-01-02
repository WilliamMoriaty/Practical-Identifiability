function dxdthe3 = partial_theta3(t,x,para)
    % \partial x/\partial \theta_1
    %[~, Y] = ode15s(@(t,y)Michaelis_Menten_eq(t,y,para), tspan1, y0);
    %y1_t = Y(end,1);y2_t = Y(end,2);y3_t = Y(end,3); y4_t = Y(end,4);
    %J = Jacobian(para,y1_t,y2_t);
    % ODE
     k1 = para(1);
    k2 = para(2);
    k3 = para(3);

    y1 = x(5);
    y2 = x(6);
    y3 = x(7);
    dxdthe3 = zeros(8,1);
    J = Jacobian(para,y1,y2);
    dxdthe3(1) = 0 + J(1,1)*x(1) + J(1,2)*x(2) + J(1,3)*x(3) + J(1,4)*x(4);
    dxdthe3(2) = y3 + J(2,1)*x(1) + J(2,2)*x(2) + J(2,3)*x(3) + J(2,4)*x(4);
    dxdthe3(3) = -y3 + J(3,1)*x(1) + J(3,2)*x(2) + J(3,3)*x(3) + J(3,4)*x(4);
    dxdthe3(4) = y3 + J(4,1)*x(1) + J(4,2)*x(2) + J(4,3)*x(3) + J(4,4)*x(4);
    dxdthe3(5) = -k1 * y1 * y2 + k2 * y3;  % dy1/dt
    dxdthe3(6) = -k1 * y1 * y2 + (k2 + k3) * y3; % dy2/dt
    dxdthe3(7) = k1 * y1 * y2 - (k2 + k3) * y3; % dy3/dt
    dxdthe3(8) = k3 * y3; % dy4/dt
end