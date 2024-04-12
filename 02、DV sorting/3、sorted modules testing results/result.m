%%  看看效果
% 同组电池的短期一致性（7天的电压和容量）和长期一致性（半年的电压和容量）；
% 与随机分配的对比。
clear;clc
cd("C:\Users\11566\Desktop\毕业论文书写\第三章\materials_new\效果\");
%%  计算分组的电压和容量的标准差
%加载数据
% 找最符合的，所有CAP和VOLT都引入。。
VOLT=load("VOLT.mat");
VOLT=VOLT.VOLT;
CAP=load("CAP.mat");
CAP=CAP.CAP;
%%  
thegroup=load('group_dv1.mat');
thegroup=thegroup.group;

for i=1:length(thegroup)
    onegroup=thegroup{i};
    %电压
    vol=VOLT(onegroup,[2 end]);
    %容量
    cap=CAP(onegroup,[2 end]);
%     if i==1
%     z1=min(cap(:,1))./mean(cap(:,1))-1;
%     z2=max(cap(:,1))./mean(cap(:,1))-1;
%     end
    ALL_CK=[vol cap];
%     ALL_CK=normalize(ALL_CK,'range');
    ALL_talk_CK(i,:)=std(ALL_CK);
    ALL_talk_mean(i,:)=mean(ALL_CK);
    jiaquan(i)=length(onegroup);
end
   jiaquan1=jiaquan./sum(jiaquan);
   sum_in_CK=jiaquan1*ALL_talk_CK;
%%
Need=[ALL_talk_CK;sum_in_CK];
% Need(:,4)=1.1*Need(:,4);
% Need(:,3)=1.1*Need(:,3);