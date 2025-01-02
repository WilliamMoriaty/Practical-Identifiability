function dydt = lotka_volterra_eq(t,y,para)

% Lotka-Volterra equations

alpha = para(1);
beta = para(2);
delta = para(3);
gamma = para(4);
    x = y(1);
    y = y(2);
    dydt = zeros(2,1);
    dydt(1) = alpha * x - beta * x * y;  % Prey equation
    dydt(2) = delta * x * y - gamma * y; % Predator equation
end