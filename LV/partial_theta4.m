
function dxdthe4 = partial_theta4(t, x,para)
    % \partial x/\partial theta_4
%     [~, Y] = ode45(@(t,y)lotka_volterra_eq(t,y,para), tspan1, y0);
%     X_t = Y(end,1);Y_t = Y(end,2);
%     J = Jacobian(para,X_t,Y_t);
    alpha = para(1);
    beta = para(2);
    delta = para(3);
    gamma = para(4);
    % ODE
    dxdthe4 = zeros(4,1);
    J = Jacobian(para,x(3),x(4));
    dxdthe4(1) = J(1,1)*x(1) + J(1,2)*x(2);
    dxdthe4(2) = -x(4) + J(2,1)*x(1) + J(2,2)*x(2);
    dxdthe4(3) = alpha * x(3) - beta * x(3) * x(4);  % Prey equation
    dxdthe4(4) = delta * x(3) * x(4) - gamma * x(4); % Predator equation

end