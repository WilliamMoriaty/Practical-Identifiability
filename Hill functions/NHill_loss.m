function L=NHill_loss(t,x1,x2,x3,y,lambda,U)
L=0;
for i=1:length(t)
L1 = (x1*(t(i).^x3./(t(i).^x3+x2.^x3))-y(i)).^2;
L = L+L1;
end
xs=[1 3 8];
x=[x1 x2 x3];
L=L+lambda*norm(U'*(x-xs)',2)^2;
end