function dxdthe2 = partial_theta2(t, x,para)
    % \partial x/\partial \theta_2
%     [~, Y] = ode45(@(t,y)lotka_volterra_eq(t,y,para), tspan1, y0);
%     X_t = Y(end,1);Y_t = Y(end,2);
%    J = Jacobian(para,X_t,Y_t);
    alpha = para(1);
    beta = para(2);
    delta = para(3);
    gamma = para(4);
%     x(3) = y(1); prey
%     x(4) = y(2); predator
%     dydt = zeros(2,1);
%     dydt(1) = alpha * x - beta * x * y;  % Prey equation
%     dydt(2) = delta * x * y - gamma * y; % Predator equation

    % ODE
    dxdthe2 = zeros(4,1);
    J = Jacobian(para,x(3),x(4));
    dxdthe2(1) = -x(3)*x(4) + J(1,1)*x(1) + J(1,2)*x(2);
    dxdthe2(2) = J(2,1)*x(1) + J(2,2)*x(2);
    dxdthe2(3) = alpha * x(3) - beta * x(3) * x(4);  % Prey equation
    dxdthe2(4) = delta * x(3) * x(4) - gamma * x(4); % Predator equation

end