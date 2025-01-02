function L=Hill_loss(t,x1,x2,x3,y)
L=0;
for i=1:length(t)
L1 = (x1*(t(i).^x3./(t(i).^x3+x2.^x3))-y(i)).^2;
L = L+L1;
end
end