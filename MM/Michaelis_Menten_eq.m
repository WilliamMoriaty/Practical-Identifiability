function dydt = Michaelis_Menten_eq(t,y,para)

%  Michaelis-Menten equations

k1 = para(1);
k2 = para(2);
k3 = para(3);

    y1 = y(1);
    y2 = y(2);
    y3 = y(3);
    
    dydt = zeros(4,1);
    dydt(1) = -k1 * y1 * y2 + k2 * y3;  % dy1/dt
    dydt(2) = -k1 * y1 * y2 + (k2 + k3) * y3; % dy2/dt
    dydt(3) = k1 * y1 * y2 - (k2 + k3) * y3; % dy3/dt
    dydt(4) = k3 * y3; % dy4/dt
end