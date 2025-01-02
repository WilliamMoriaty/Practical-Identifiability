function L=NLinear_loss(t,x1,x2,x3,y,lambda,U)
L=0;
for i=1:length(t)
L1 = (x1+x2*t(i)^2+x3*((t(i)-1)*(t(i)-2)*(t(i)-3)+2)-y(i)).^2;
L = L+L1;
end
xs=[1 1 1];
x=[x1 x2 x3];
L=L+lambda*norm(U'*(x-xs)',2)^2;
end