function L=Linear_loss(t,x1,x2,x3,y)
L=0;
for i=1:length(t)
L1 = (x1+x2*t(i)^2+x3*((t(i)-1)*(t(i)-2)*(t(i)-3)+2)-y(i)).^2;
L = L+L1;
end
end