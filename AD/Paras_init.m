% dt=1e-1;
% ParasValues_AD
% opts=odeset('NonNegative',1,'AbsTol',1e-8,'RelTol',1e-8);% Impose nonnegative constraint on parameters
%% AD Data discovery
AD_ind = find(txt{2}(:,8)=="AD");
% LMCI_ind = find(txt{2}(:,8)=="LMCI");
%  CN_ind = find(txt{2}(:,8)=="CN");

flag_mat=[];Ndata=[];
for j=1:length(AD_ind)
flag=0;
% ind=2;
% PID=find(data{4}(:,1)==patient_id(ind));
PID=find(data{4}(:,1)==AD_ind(j));
Abeta_PID=Abeta(PID);
ptau_PID=ptau(PID);
tau0_PID=tau0(PID);
ref_date = datetime(1900, 1, 1);  % date baseline
AAge=zeros(size(PID));
select_id=[];
for i=1:length(PID)
    %AAge(i)=Age(ind)+(datenum(txt{4}(PID(i)+1,6))-Age_time(ind))/365;
    ind=find(patient_id==AD_ind(j));
    AAge(i)=Age(ind)+((data{4}(PID(i),6))-days(Age_time(ind) - ref_date))/365;
    if strcmp(txt{4}(PID(i)+1,3),'MEDIAN')
        select_id=[select_id i];
    end
end
Abeta_PID=Abeta_PID(select_id);
ptau_PID=ptau_PID(select_id);
tau0_PID=tau0_PID(select_id);

AAge=AAge(select_id);
pAge=AAge;
tAge=AAge;
AAge(isnan(Abeta_PID))=[];
Abeta_PID(isnan(Abeta_PID))=[];

pAge(isnan(ptau_PID))=[];
ptau_PID(isnan(ptau_PID))=[];

tAge(isnan(tau0_PID))=[];
tau0_PID(isnan(tau0_PID))=[];

if length(Abeta_PID) <4 || length(tau0_PID) <4 || length(ptau_PID) <4
    flag=-1;
%     return;
else
% xxx = [length(Abeta_PID) length(tau0_PID) length(ptau_PID)];
%  pathname = strcat("AD_",num2str(AD_ind(j)),'.mat');
%  save(pathname,'Abeta_PID','tau0_PID','ptau_PID','AAge')
%  pathname = strcat("LMCI_",num2str(LMCI_ind(j)),'.mat');
%  save(pathname,'Abeta_PID','tau0_PID','ptau_PID','AAge')
pathname = strcat("AD_data/AD_",num2str(AD_ind(j)),'_13.mat');
save(pathname,'Abeta_PID','tau0_PID','ptau_PID','AAge')
%Ndata = [Ndata;xxx];
end
% flag_mat = [flag_mat flag]; 
end

for j=1:length(AD_ind)
        PID=find(data{4}(:,1)==AD_ind(j));
        N_PID=N(PID);
        C_PID=C(PID);
% 
        NAge=zeros(size(PID));
        select_id=[];

        for i=1:length(PID)
            %NAge(i)=Age(ind)+(datenum(txt{2}(PID(i)+1,7))-Age_time(ind))/365;
            ind=find(patient_id==AD_ind(j));
            NAge(i)=Age(ind)+((data{4}(PID(i),7))-days(Age_time(ind) - ref_date))/365;
    
            if strcmp(txt{2}(PID(i)+1,5),'ADNI1')
                select_id=[select_id i];
            end
        end
        N_PID=N_PID(select_id);
        C_PID=C_PID(select_id);
        NAge=NAge(select_id);
        CAge=NAge;
        CAge(isnan(C_PID))=[];
        C_PID(isnan(C_PID))=[];
        NAge(isnan(N_PID))=[];
        N_PID(isnan(N_PID))=[];
        if length(N_PID) <3 || length(C_PID) <3
            flag=-1;
        else
        
        pathname = strcat("AD_data/AD_",num2str(AD_ind(j)),'_45.mat');
        save(pathname,'N_PID','C_PID','NAge')
        end
        
end
%% LMCI Data discovery
%AD_ind = find(txt{2}(:,8)=="AD");
 LMCI_ind = find(txt{2}(:,8)=="LMCI");
%  CN_ind = find(txt{2}(:,8)=="CN");

flag_mat=[];Ndata=[];
for j=1:length(LMCI_ind)
flag=0;
% ind=2;
% PID=find(data{4}(:,1)==patient_id(ind));
PID=find(data{4}(:,1)==LMCI_ind(j));
Abeta_PID=Abeta(PID);
ptau_PID=ptau(PID);
tau0_PID=tau0(PID);
ref_date = datetime(1900, 1, 1);  % date baseline
AAge=zeros(size(PID));
select_id=[];
for i=1:length(PID)
    %AAge(i)=Age(ind)+(datenum(txt{4}(PID(i)+1,6))-Age_time(ind))/365;
    ind=find(patient_id==LMCI_ind(j));
    AAge(i)=Age(ind)+((data{4}(PID(i),6))-days(Age_time(ind) - ref_date))/365;
    if strcmp(txt{4}(PID(i)+1,3),'MEDIAN')
        select_id=[select_id i];
    end
end
Abeta_PID=Abeta_PID(select_id);
ptau_PID=ptau_PID(select_id);
tau0_PID=tau0_PID(select_id);

AAge=AAge(select_id);
pAge=AAge;
tAge=AAge;
AAge(isnan(Abeta_PID))=[];
Abeta_PID(isnan(Abeta_PID))=[];

pAge(isnan(ptau_PID))=[];
ptau_PID(isnan(ptau_PID))=[];

tAge(isnan(tau0_PID))=[];
tau0_PID(isnan(tau0_PID))=[];

if length(Abeta_PID) <6 || length(tau0_PID) <6 || length(ptau_PID) <6
    flag=-1;
%     return;
else
% xxx = [length(Abeta_PID) length(tau0_PID) length(ptau_PID)];
%  pathname = strcat("AD_",num2str(AD_ind(j)),'.mat');
%  save(pathname,'Abeta_PID','tau0_PID','ptau_PID','AAge')
%  pathname = strcat("LMCI_",num2str(LMCI_ind(j)),'.mat');
%  save(pathname,'Abeta_PID','tau0_PID','ptau_PID','AAge')
pathname = strcat("LMCI_data/LMCI_",num2str(LMCI_ind(j)),'_13.mat');
save(pathname,'Abeta_PID','tau0_PID','ptau_PID','AAge')
%Ndata = [Ndata;xxx];
end
% flag_mat = [flag_mat flag]; 
end

for j=1:length(LMCI_ind)
        PID=find(data{4}(:,1)==LMCI_ind(j));
        N_PID=N(PID);
        C_PID=C(PID);
% 
        NAge=zeros(size(PID));
        select_id=[];

        for i=1:length(PID)
            %NAge(i)=Age(ind)+(datenum(txt{2}(PID(i)+1,7))-Age_time(ind))/365;
            ind=find(patient_id==LMCI_ind(j));
            NAge(i)=Age(ind)+((data{4}(PID(i),7))-days(Age_time(ind) - ref_date))/365;
    
            if strcmp(txt{2}(PID(i)+1,5),'ADNI1')
                select_id=[select_id i];
            end
        end
        N_PID=N_PID(select_id);
        C_PID=C_PID(select_id);
        NAge=NAge(select_id);
        CAge=NAge;
        CAge(isnan(C_PID))=[];
        C_PID(isnan(C_PID))=[];
        NAge(isnan(N_PID))=[];
        N_PID(isnan(N_PID))=[];
        if length(N_PID) <3 || length(C_PID) <3
            flag=-1;
        else
       
        pathname = strcat("LMCI_data/LMCI_",num2str(LMCI_ind(j)),'_45.mat');
        save(pathname,'N_PID','C_PID','NAge')
        end
        
end


%% CN Data discovery
%AD_ind = find(txt{2}(:,8)=="AD");
% LMCI_ind = find(txt{2}(:,8)=="LMCI");
  CN_ind = find(txt{2}(:,8)=="CN");

flag_mat=[];Ndata=[];
for j=1:length(CN_ind)
flag=0;
% ind=2;
% PID=find(data{4}(:,1)==patient_id(ind));
PID=find(data{4}(:,1)==CN_ind(j));
Abeta_PID=Abeta(PID);
ptau_PID=ptau(PID);
tau0_PID=tau0(PID);
ref_date = datetime(1900, 1, 1);  % date baseline
AAge=zeros(size(PID));
select_id=[];
for i=1:length(PID)
    %AAge(i)=Age(ind)+(datenum(txt{4}(PID(i)+1,6))-Age_time(ind))/365;
    ind=find(patient_id==CN_ind(j));
    AAge(i)=Age(ind)+((data{4}(PID(i),6))-days(Age_time(ind) - ref_date))/365;
    if strcmp(txt{4}(PID(i)+1,3),'MEDIAN')
        select_id=[select_id i];
    end
end
Abeta_PID=Abeta_PID(select_id);
ptau_PID=ptau_PID(select_id);
tau0_PID=tau0_PID(select_id);

AAge=AAge(select_id);
pAge=AAge;
tAge=AAge;
AAge(isnan(Abeta_PID))=[];
Abeta_PID(isnan(Abeta_PID))=[];

pAge(isnan(ptau_PID))=[];
ptau_PID(isnan(ptau_PID))=[];

tAge(isnan(tau0_PID))=[];
tau0_PID(isnan(tau0_PID))=[];

if length(Abeta_PID) <6 || length(tau0_PID) <6 || length(ptau_PID) <6
    flag=-1;
%     return;
else
% xxx = [length(Abeta_PID) length(tau0_PID) length(ptau_PID)];
%  pathname = strcat("AD_",num2str(AD_ind(j)),'.mat');
%  save(pathname,'Abeta_PID','tau0_PID','ptau_PID','AAge')
%  pathname = strcat("LMCI_",num2str(LMCI_ind(j)),'.mat');
%  save(pathname,'Abeta_PID','tau0_PID','ptau_PID','AAge')
pathname = strcat("CN_data/CN_",num2str(CN_ind(j)),'_13.mat');
save(pathname,'Abeta_PID','tau0_PID','ptau_PID','AAge')
%Ndata = [Ndata;xxx];
end
% flag_mat = [flag_mat flag]; 
end


for j=1:length(CN_ind)
        PID=find(data{4}(:,1)==CN_ind(j));
        N_PID=N(PID);
        C_PID=C(PID);
% 
        NAge=zeros(size(PID));
        select_id=[];

        for i=1:length(PID)
            %NAge(i)=Age(ind)+(datenum(txt{2}(PID(i)+1,7))-Age_time(ind))/365;
            ind=find(patient_id==CN_ind(j));
            NAge(i)=Age(ind)+((data{4}(PID(i),7))-days(Age_time(ind) - ref_date))/365;
    
            if strcmp(txt{2}(PID(i)+1,5),'ADNI1')
                select_id=[select_id i];
            end
        end
        N_PID=N_PID(select_id);
        C_PID=C_PID(select_id);
        NAge=NAge(select_id);
        CAge=NAge;
        CAge(isnan(C_PID))=[];
        C_PID(isnan(C_PID))=[];
        NAge(isnan(N_PID))=[];
        N_PID(isnan(N_PID))=[];
        if length(N_PID) <3 || length(C_PID) <3
            flag=-1;
        else

        pathname = strcat("CN_data/CN_",num2str(CN_ind(j)),'_45.mat');
        save(pathname,'N_PID','C_PID','NAge')
        end       
end
%%
% %%
% a=min(Abeta_PID);
% %60
% x=fmincon(@(x) Parameter_est_obj_A(x,Abeta_PID,AAge,paras),[0.1;300;10],[],[],[],[],[0;0;0],[0.2;400;60]);
% u0=paras.initial;
% paras.lambda_beta=x(1);
% paras.K_Abeta=x(2);
% u0(1)=x(3);
% paras.initial=u0;
% [t,u]=ode45(@(t,u) Right_hand_side(t,u,paras),[50:0.1:100],u0,opts);
% n=length(Abeta_PID);
% f_A=sqrt(1/n*sum(((u(floor((AAge-50)/0.1),1)-Abeta_PID)./Abeta_PID).^2));
% if f_A<Tol
%     a=min(ptau_PID);
%     x=fmincon(@(x) Parameter_est_obj_p(x,ptau_PID,pAge,paras),[0.02;60;0.1],[],[],[],[],[0;0;0],[0.5;200;a]);
%     paras.lambda_ptau=x(1);
%     paras.K_ptau=x(2);
%     u0(2)=x(3);
%     paras.initial=u0;
% 
%     [t,u]=ode45(@(t,u) Right_hand_side(t,u,paras),[50:0.1:100],u0,opts);
%     n=length(ptau_PID);
%     f_p=sqrt(1/n*sum(((u(floor((pAge-50)/0.1),2)-ptau_PID)./ptau_PID).^2));
%     if f_p<Tol
%         a=min(tau0_PID);
%         x=fmincon(@(x) Parameter_est_obj_t(x,tau0_PID,tAge,paras),[0;20],[],[],[],[],[0;0],[10;a]);
%         paras.lambdatauO=x(1);
%         u0(3)=x(2);
%         paras.initial=u0;
%         [t,u]=ode45(@(t,u) Right_hand_side(t,u,paras),[50:0.1:100],u0,opts);
%         n=length(tau0_PID);
%         f_o=sqrt(1/n*sum(((u(floor((tAge-50)/0.1),3)-tau0_PID)./tau0_PID).^2));
% 
%         PID=find(data{2}(:,1)==patient_id(ind));
%         N_PID=N(PID);
%         C_PID=C(PID);
% % 
% %         NAge=zeros(size(PID));
% %         select_id=[];
% % 
% %         for i=1:length(PID)
% %             NAge(i)=Age(ind)+(datenum(txt{2}(PID(i)+1,7))-Age_time(ind))/365;
% %             if strcmp(txt{2}(PID(i)+1,5),'ADNI1')
% %                 select_id=[select_id i];
% %             end
% %         end
% %         N_PID=N_PID(select_id);
% %         C_PID=C_PID(select_id);
% %         NAge=NAge(select_id);
% %         CAge=NAge;
% %         CAge(isnan(C_PID))=[];
% %         C_PID(isnan(C_PID))=[];
% %         NAge(isnan(N_PID))=[];
% %         N_PID(isnan(N_PID))=[];
% %         if length(N_PID) <3 || length(C_PID) <3
% %             flag=-1;
% %             return;
% %         end
% % 
%         opts1 = optimoptions('fmincon','StepTolerance',1e-30,'OptimalityTolerance',1e-30);
%         x=fmincon(@(x) Parameter_est_obj_N(x,N_PID,NAge,paras),[0.0001;0.01;1;0.001],[],[],[],[],[0;0;0;0],[0.0005;0.008;1.2;0.3],[],opts1);
%         paras.lambda_NdTO=x(1);
%         paras.lambda_Ndptau=x(2);
%         paras.K_Nd=x(3);
%         u0(4)=x(4);
%         paras.initial=u0;
%         [t,u]=ode45(@(t,u) Right_hand_side(t,u,paras),[50:0.1:100],u0,opts);
%         n=length(N_PID);
%         f_N=sqrt(1/n*sum(((u(floor((NAge-50)/0.1),4)-N_PID)./N_PID).^2));
% 
%         if f_N<Tol
%             a=min(C_PID);
%             %15
%             b=max(C_PID);
%             %50
%             x=fmincon(@(x) Parameter_est_obj_C(x,C_PID,CAge,paras),[0.05;0.1;50;1],[],[],[],[],[0;0;0;0],[10;20;200;a]);
%             paras.lambda_CN=x(1);
%             paras.lambda_Ct=x(2);
%             paras.K_C=x(3);
%             u0(5)=x(4);
%             paras.initial=u0;
% 
%             [t,u]=ode45(@(t,u) Right_hand_side(t,u,paras),[50:0.1:100],u0,opts);
%             n=length(C_PID);
%             f_C=sqrt(1/n*sum(((u(floor((CAge-50)/0.1),5)-C_PID)./C_PID).^2));
% 
%             if f_C<Tol
% 
%                 x0(1)=paras.lambda_beta;
%                 x0(2)=paras.K_Abeta;
%                 x0(3)=u0(1);
% 
%                 x0(4)=paras.lambda_ptau;
%                 x0(5)=paras.K_ptau;
%                 x0(6)=u0(2);
% 
%                 x0(7)=paras.lambdatauO;
%                 x0(8)=u0(3);
% 
%                 x0(9)=paras.lambda_NdTO;
%                 x0(10)=paras.lambda_Ndptau;
%                 x0(11)=paras.K_Nd;
%                 x0(12)=u0(4);
% 
%                 x0(13)=paras.lambda_CN;
%                 x0(14)=paras.lambda_Ct;
%                 x0(15)=paras.K_C;
%                 x0(16)=u0(5);
%             else
%                 flag=-2;
%                 return;
% 
%             end
%         else
%             flag=-2;
%             return;
%         end
%     else
%         flag=-2;
%         return;
%     end
% else 
%     flag=-2;
%     return;
% end
