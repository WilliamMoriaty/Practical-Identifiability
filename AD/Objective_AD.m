function f = Objective_AD(para,tspan1,tspan2,y0,yexp1,yexp2,yexp3,yexp4,lambda,U,para0)
       f=0;
        sol = ode45(@(t,y)AD_eq(t,y,para),tspan1,y0);
        solpts = deval(sol,tspan1);
        f = f+sum((solpts(1,:) - yexp1').^2)/size(yexp1,1); %  y(3)
    
        sol = ode45(@(t,y)AD_eq(t,y,para),tspan1,y0);
        solpts = deval(sol,tspan1);
        f = f + sum((solpts(2,:) - yexp2').^2)/size(yexp2,1);
    
        sol = ode45(@(t,y)AD_eq(t,y,para),tspan2,y0);
        solpts = deval(sol,tspan2);
        f = f + sum((solpts(3,:) - yexp3').^2)/size(yexp3,1);
%   
        sol = ode45(@(t,y)AD_eq(t,y,para),tspan2,y0);
        solpts = deval(sol,tspan2);
        f = f + sum((solpts(4,:) - yexp4').^2)/size(yexp4,1);


        f = f + lambda*norm(U'*(para-para0),2)*norm(U'*(para-para0),2);
end
