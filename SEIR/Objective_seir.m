function f = Objective_seir(para,tspan,y0,index,yexp1,yexp2,lambda,U,para0)

switch index
    case 1
        sol = ode45(@(t,y)SEIR_eq(t,y,para),tspan,y0);
        solpts = deval(sol,tspan);
        f = sum((solpts(3,:) - yexp1').^2)/size(yexp1,1); %  y(3)
    case 2
        sol = ode45(@(t,y)SEIR_eq(t,y,para),tspan,y0);
        solpts = deval(sol,tspan);
        f = sum((solpts(4,:) - yexp2').^2)/size(yexp2,1);
    case 3
        sol = ode45(@(t,y)SEIR_eq(t,y,para),tspan,y0);
        solpts = deval(sol,tspan);
        f = sum((solpts(3,:) - yexp1').^2)/size(yexp1,1);
        f = f+sum((solpts(4,:) - yexp2').^2)/size(yexp2,1);
    case 4
        sol = ode45(@(t,y)SEIR_eq(t,y,para),tspan,y0);
        solpts = deval(sol,tspan);
        f = sum((solpts(3,:) - yexp1').^2)/size(yexp1,1); %  y(3)
        f = f + lambda*norm(U'*(para-para0),2)*norm(U'*(para-para0),2);
end

%ff=f1+f2+f3+lambda*norm(par,2);
end