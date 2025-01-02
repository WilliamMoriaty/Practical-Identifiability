[~,sheet_name]=xlsfinfo('ADNIMERGE Customized CSF.xlsx');
for k=1:numel(sheet_name)
    [data{k},txt{k}]=xlsread('ADNIMERGE Customized CSF.xlsx',sheet_name{k});
end
[patient_id,ind]=unique(data{2}(:,1));
Time= data{2}(:,7);
%date_times = datetime(Time, 'ConvertFrom', 'excel');

Age= data{2}(ind,9);
Gender=txt{2}(ind+1,10);
%formatIn = 'yyyy/mm/dd';
Age_time=Age;
% Initialize Age_time 
%Age_time = NaT(length(patient_id), 1);  % NaT 表示 "Not-a-Time"

for i=1:length(patient_id)
    ind=find(data{2}(:,1)==patient_id(i));
    Age_time(i)=min(data{2}(ind,7));
    %Age_time(i)=min(datenum(txt{2}(ind+1,7)));
    %Age_time(i) = min(datetime(txt{2}(ind + 1, 7), 'InputFormat', 'yyyy/MM/dd'));  % 使用 datetime
end
Age_time = datetime(Age_time, 'ConvertFrom', 'excel');

Abeta= data{4}(:,8);
Abetaa=min(Abeta);
Abetab=max(Abeta);
Abeta=Abetab-Abeta;

ptau= data{4}(:,10);
ptaua=min(ptau);
ptaub=max(ptau);

tau0= data{4}(:,9);
tau0=tau0-ptau;
tau0b=max(tau0);


N= data{2}(:,46);
whole=data{2}(:,47);
N=1-N./whole;
Na=min(N);
Nb=max(N);



C= data{2}(:,21);
