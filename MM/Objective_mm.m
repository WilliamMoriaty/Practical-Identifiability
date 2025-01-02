function f = Objective_mm(para,tspan,y0,index,yexp1,yexp2)

switch index
    case 1
        sol = ode15s(@(t,y)Michaelis_Menten_eq(t,y,para),tspan,y0);
        solpts = deval(sol,tspan);
        f = sum((solpts(1,:) - yexp1').^2)/size(yexp1,1); %  y(3)
    case 2
        sol = ode15s(@(t,y)Michaelis_Menten_eq(t,y,para),tspan,y0);
        solpts = deval(sol,tspan);
        f = sum((solpts(4,:) - yexp2').^2)/size(yexp2,1);
    case 3
        sol = ode15s(@(t,y)Michaelis_Menten_eq(t,y,para),tspan,y0);
        solpts = deval(sol,tspan);
        f = sum((solpts(1,:) - yexp1').^2)/size(yexp1,1);
        f = f+sum((solpts(4,:) - yexp2').^2)/size(yexp2,1);
end

%ff=f1+f2+f3+lambda*norm(par,2);
end