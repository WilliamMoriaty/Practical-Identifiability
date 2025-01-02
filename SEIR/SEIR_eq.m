function dydt = SEIR_eq(t,y,para)

% SEIR equations

beta = para(1);
sigma = para(2);
gamma = para(3);

    S = y(1);
    E = y(2);
    I = y(3);
    
    dydt = zeros(4,1);
    dydt(1) = -beta * S * I ;  % dy1/dt
    dydt(2) = beta * S * I  - sigma * E; % dy2/dt
    dydt(3) = sigma * E - gamma * I; % dy3/dt
    dydt(4) = gamma * I; % dy4/dt
end