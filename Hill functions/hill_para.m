
function yy=hill_para(x,Vmax,Kd,n)

yy=zeros(1,3);
yy(1)=x.^n./(x.^n+Kd.^n);
yy(2)=-Vmax*n*Kd.^(n-1)*x.^n./(x.^n+Kd.^n).^2;
yy(3)=Vmax*Kd.^n*x.^n.*(log(x)-log(2))./(x.^n+Kd.^n).^2;
end