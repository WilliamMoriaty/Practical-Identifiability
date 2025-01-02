function dxdthe9 = partial_theta9(t,x,para)
    % \partial x/\partial \theta_9

    lambda_Ab = para(1);
    lambda_tau = para(2);
    lambda_Ntaup = para(3);
    lambda_CN = para(4);
    lambda_Ctau = para(5);
    K_Ab = para(6);
    K_taup = para(7);
    K_N = para(8);
    K_C = para(9);

    Ab = x(5);
    tau_p = x(6);
    N = x(7);
    C = x(8);
    
    % ODE
    dxdthe9 = zeros(8,1);
    J = Jacobian(para,Ab,tau_p,N,C);
    dxdthe9(1) = 0 + J(1,1)*x(1) + J(1,2)*x(2) + J(1,3)*x(3) + J(1,4)*x(4);
    dxdthe9(2) = 0 + J(2,1)*x(1) + J(2,2)*x(2) + J(2,3)*x(3) + J(2,4)*x(4);
    dxdthe9(3) = 0 + J(3,1)*x(1) + J(3,2)*x(2) + J(3,3)*x(3) + J(3,4)*x(4);
    dxdthe9(4) = (lambda_CN*N+lambda_Ctau*tau_p)*C/(K_C*K_C) + J(4,1)*x(1) + J(4,2)*x(2) + J(4,3)*x(3) + J(4,4)*x(4);
    dxdthe9(5) = lambda_Ab * Ab * (1 - Ab/K_Ab) ;  % dAb/dt
    dxdthe9(6) = lambda_tau * Ab * (1 - tau_p/K_taup); % dtau_p/dt
    dxdthe9(7) = (lambda_Ntaup * tau_p) * (1 - N/K_N); % dN/dt
    dxdthe9(8) = (lambda_CN * N + lambda_Ctau * tau_p) * (1 - C/K_C); % dC/dt

end