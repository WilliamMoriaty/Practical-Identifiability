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
for m=m1+1:M+m1
T = max(t)+dt;

yy=hill_para(T,Vmax,Kd,n);
x = diag((yy*U(:,end))'*yy*U(:,end));
non_zero_indices = x ~= 0;
num_non_zero = sum(non_zero_indices);

if num_non_zero > 0
    break;
end


end