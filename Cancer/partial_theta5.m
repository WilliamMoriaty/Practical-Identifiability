function dxdthe5 = partial_theta5(t,x,para,n)
    % \partial x/\partial \theta_5
    dx=1/n;D1 = 1e-6;D2=1e-6;
    sigma = para(1);rho = para(2);eta = para(3);
    mu = para(4);eps = para(5);
    beta1 = para(6);beta2 = para(7);
    phi = para(8);lambda = para(9);
    psi = para(10);omega=D2/D1;
    E = x(3*(n+1)+1:4*(n+1));
    T = x(4*(n+1)+1:5*(n+1));
    C = x(5*(n+1)+1:6*(n+1));
    % ODE
    dxdthe5 = zeros(3*(n+1)*2,1);
    for i = 1:n+1
    J = Jacobian(para,E(i),T(i),C(i));
    ff1 = C(i) + J(1,1)*x(i) + J(1,2)*x(n+1+i) + J(1,3)*x(2*(n+1)+i);
    ff2 = 0 + J(2,1)*x(i) + J(2,2)*x(n+1+i) + J(2,3)*x(2*(n+1)+i);
    ff3 = 0 + J(3,1)*x(i) + J(3,2)*x(n+1+i) + J(3,3)*x(2*(n+1)+i);
    if i == 1
    dxdthe5(i) = (2/dx^2)*(x(i+1)-x(i))+ ff1;
    dxdthe5(n+1+i) = (2*omega/dx^2)*(x(n+1+1+i)-x(n+1+i)) + ff2;
    dxdthe5(2*(n+1)+i) = ff3;
    elseif i == n+1
    dxdthe5(i) = (2/dx^2)*(x(i-1)-x(i)) + ff1;
    dxdthe5(n+1+i) = (2*omega/dx^2)*(x(n+i)-x(n+1+i)) + ff2;
    dxdthe5(2*(n+1)+i) = ff3;
    else
    dxdthe5(i) = (1/dx^2)*(x(i+1)+x(i-1)-2*x(i)) + ff1;
    dxdthe5(n+1+i) = (omega/dx^2)*(x(n+i+1+1)+x(n+i)-2*x(n+i+1)) + ff2;
    dxdthe5(2*(n+1)+i) = ff3;
    end

    end

    for i = 1:n+1
        fE = sigma*hvd((i-1)*dx)+rho*C(i)/(eta+T(i))-sigma*E(i)-mu*E(i)*T(i)+eps*C(i);
        fT = beta1*T(i)*(1-beta2*T(i))-phi*E(i)*T(i)+lambda*C(i);
        fC = mu*E(i)*T(i)-psi*C(i);
       
    if i == 1
       dxdthe5(3*(n+1)+i)=(2/dx^2)*(E(i+1)-E(i))+fE;
       dxdthe5(4*(n+1)+i)=(2*omega/dx^2)*(T(i+1)-T(i))+fT;
       dxdthe5(5*(n+1)+i)=fC;
    elseif i == n+1
       dxdthe5(3*(n+1)+i)=(2/dx^2)*(E(i-1)-E(i))+fE;
       dxdthe5(4*(n+1)+i)=(2*omega/dx^2)*(T(i-1)-T(i))+fT;
       dxdthe5(5*(n+1)+i)=fC;
    else
       dxdthe5(3*(n+1)+i)=(1/dx^2)*(E(i+1)+E(i-1)-2*E(i))+fE;
       dxdthe5(4*(n+1)+i)=(omega/dx^2)*(T(i+1)+T(i-1)-2*T(i))+fT;
       dxdthe5(5*(n+1)+i)=fC;
    end

   end

end
