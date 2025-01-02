t = [1.6,3.2]; 
Vmax = 1.0;  % 
Kd = 3.0;    % 
n = 8;       % 

y = hill_function(t, Vmax, Kd, n);

yy=zeros(size(t,2),3);
for i=1:size(t,2)
yy(i,:)=hill_para(t(i),Vmax,Kd,n);
end
yy3=yy;
[U,Sigma,~]=svd(yy3'*yy3);

dt=2.6;
eps = 1e-6;
M = 5;
m1 = size(t,2);
alpha = zeros(1,3);
sigmak = Sigma(3,3);
for m=m1:M+m1
T = max(t)+dt;
t = [t,T];

for i=1:size(t,2)
yy(i,:)=hill_para(t(i),Vmax,Kd,n);
end

for j=1:3
alpha(j)=yy(j,:)*U(:,j);
end
if size(find(alpha==0)) > 0
    continue;
end

[~,Sigma,~]=svd(yy'*yy);
sigmak = Sigma(3,3);
if sigmak > eps
    break;
end

end