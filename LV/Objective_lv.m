function f = Objective_lv(para,tspan,y0,index,yexp1,yexp2,lambda,U,para0)

switch index
    case 1
        sol = ode45(@(t,y)lotka_volterra_eq(t,y,para),tspan,y0);
        solpts = deval(sol,tspan);
        f = sum((solpts(1,:) - yexp1').^2)/size(yexp1,1);
        f = f+sum((solpts(2,:) - yexp2').^2)/size(yexp2,1);% y(:,1:2)
    case 2
        %lambda = 0.1;
        sol = ode45(@(t,y)lotka_volterra_eq(t,y,para),tspan,y0);
        solpts = deval(sol,tspan);
        f = sum((solpts(1,:) - yexp1').^2)/size(yexp1,1);
        f = f+sum((solpts(2,:) - yexp2').^2)/size(yexp2,1);% y(:,1:2)
        f = f + lambda*norm(U'*(para-para0),2)*norm(U'*(para-para0),2);
end

%ff=f1+f2+f3+lambda*norm(par,2);
end