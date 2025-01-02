function dydt = PDE_eq(t,y,para,n)

% PDE equations

omega = 1;
dx = 1/n;
sigma = para(1);rho = para(2);eta = para(3);
mu = para(4);eps = para(5);
beta1 = para(6);beta2 = para(7);
phi = para(8);lambda = para(9);
psi = para(10);


    E = y(1:n+1);
    T = y(n+2:2*(n+1));
    C = y(2*n+3:3*(n+1));
    
    dydt = zeros(3*(n+1),1);
   for i = 1:n+1
        fE = sigma*hvd((i-1)*dx)+rho*C(i)/(eta+T(i))-sigma*E(i)-mu*E(i)*T(i)+eps*C(i);
        fT = beta1*T(i)*(1-beta2*T(i))-phi*E(i)*T(i)+lambda*C(i);
        fC = mu*E(i)*T(i)-psi*C(i);
        

    if i == 1
       dydt(i)=(2/dx^2)*(E(i+1)-E(i))+fE;
       dydt(n+1+i)=(2*omega/dx^2)*(T(i+1)-T(i))+fT;
       dydt(2*(n+1)+i)=fC;
    elseif i == n+1
       dydt(i)=(2/dx^2)*(E(i-1)-E(i))+fE;
       dydt(n+1+i)=(2*omega/dx^2)*(T(i-1)-T(i))+fT;
       dydt(2*(n+1)+i)=fC;
    else
       dydt(i)=(1/dx^2)*(E(i+1)+E(i-1)-2*E(i))+fE;
       dydt(n+1+i)=(omega/dx^2)*(T(i+1)+T(i-1)-2*T(i))+fT;
       dydt(2*(n+1)+i)=fC;
    end

   end
end

function h=hvd(x)
if x < 0.2
   h = 0;
else
   h = 1;
end
end