Read_data
paras_list=[];
Tol=.2;
for ind=1:length(patient_id)
    Paras_init
    if flag==-1 || flag==-2
        continue
    end
    a=min(Abeta_PID);
    b=min(ptau_PID);
    c=min(tau0_PID);
    d=min(N_PID);
    e=min(C_PID);


    lb=zeros(16,1);
    ub=inf(16,1);
    ub([3,6,8,12,16])=[a;b;c;d;e];
    ub=    [0.2;400;a;0.5;200;b;10;c;0.0005;0.008;1.2;d;10;20;200;e];
    opts=optimoptions('fmincon','MaxFunctionEvaluations',100,'MaxIterations',10);
    x=fmincon(@(x) Parameter_est_obj(x,Abeta_PID,ptau_PID,tau0_PID,N_PID,C_PID,AAge,pAge,tAge,NAge,CAge,0),...
        x0,[],[],[],[],lb,ub,[],opts);

    [f,RE]=Parameter_est_obj(x,Abeta_PID,ptau_PID,tau0_PID,N_PID,C_PID,AAge,pAge,tAge,NAge,CAge,0);

    if strcmp(txt{2}(PID(1)+1,8),'AD')==1
        disease=1;
    end
    if strcmp(txt{2}(PID(1)+1,8),'CN')==1
        disease=2;
    end
    if strcmp(txt{2}(PID(1)+1,8),'LMCI')==1
        disease=3;
    end
    paras_list=[paras_list [x';RE';f;ind;disease]];
end

for draw=1:3
    disease_id=find(paras_list(end,:)==draw);
    [~,id]=min(paras_list(end-2,disease_id));
    ind=paras_list(end-1,disease_id(id));
    patient_id(ind)
    Paras_init
    x=paras_list(1:15,disease_id(id));
    [f,RE]=Parameter_est_obj(x0,Abeta_PID,ptau_PID,tau0_PID,N_PID,C_PID,AAge,pAge,tAge,NAge,CAge,draw);
end

